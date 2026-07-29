// general purpose functions for the project

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <glob.h>

#include "common.h"
#include "struct.h"
#include "auxi.h"

// fpo is a starting frequency of the band in Hz
double fpo(int band, double overlap, double dt)
{
    double fpo = 10. + (1. - overlap)*band/(2*dt);
    return fpo;
}

// fdot_min and fdot_max are the minimum and maximum frequency derivatives in
// physical units of Hz/s
void get_fdot_range(double fpo, double dt, double *fdot_min, double *fdot_max)
{
    double B = 1./(2*dt);
    // standard all-sky search range
    if (fpo < 200.) {
        *fdot_min = 2.*(fpo + B)/(2.*1000.*C_YEARSEC);
        *fdot_max = 0.;
    } else {
        *fdot_min = 2e-10;
        *fdot_max = 2e-11;
    }
}


/* Reads the input time series of detector det into ifo[det].sig.xDat,
 * counts the null samples and estimates the variance.
 *
 * The file name ifo[det].xdatname is constructed in detectors_settings()
 * (settings.c), while looking for the detector subdirectories.
 */
void read_xdat( Search_settings *sett, Command_line_opts *opts, int det )
{
     FILE *data;
     size_t status;
     int j, Nzeros=0;

     ifo[det].sig.xDat = (float *) calloc(sett->N, sizeof(float));

     if ((data = fopen(ifo[det].xdatname, "r")) != NULL) {
          if (opts->mods && strstr(opts->mods, "read_O3") != NULL) {
               // "read_O3" is present in opts->mods
               double *tmp_xdat;
               tmp_xdat = (double *) calloc(sett->N, sizeof(double));
               status = fread((void *)(tmp_xdat), sizeof(double), sett->N, data);
               for (j=0; j<sett->N; j++)
                    ifo[det].sig.xDat[j] = (float) tmp_xdat[j];
               free(tmp_xdat);
          } else {
               status = fread((void *)(ifo[det].sig.xDat), sizeof(float), sett->N, data);
          }
          fclose (data);
     } else {
          perror (ifo[det].xdatname);
          exit(EXIT_FAILURE);
     }

     // Checking for null values in the data
     for (j=0; j < sett->N; j++)
          if(!ifo[det].sig.xDat[j]) Nzeros++;

     ifo[det].sig.Nzeros = Nzeros;

     // factor N/(N - Nzeros) to account for null values in the data
     ifo[det].sig.crf0 = (double)sett->N/(sett->N - ifo[det].sig.Nzeros);

     // Estimation of the variance for each detector
     ifo[det].sig.sig2 = (ifo[det].sig.crf0)*var(ifo[det].sig.xDat, sett->N);

} // end of read xdat


/* Reads the ephemeris of detector det: its position w.r.t. the Solar System
 * barycenter for every datapoint (ifo[det].sig.DetSSB), the Earth's diurnal
 * phase phir and axis inclination epsm at t=0, and precomputes the sines and
 * cosines of the latter two.
 */
void read_detssb( Search_settings *sett, Command_line_opts *opts, int det )
{
     FILE *data;
     size_t status;
     char filename[FNAME_LENGTH];

     ifo[det].sig.DetSSB = (double *) calloc(3*sett->N, sizeof(double));
     /*
     const size_t array_bytes = 3*sett->N*sizeof(double);
     ifo[det].sig.DetSSB = NULL;
     if ( posix_memalign((void**)&ifo[det].sig.DetSSB, 32, array_bytes) ) exit (1);
     */

     sprintf (filename, "%s/%03d/%s/DetSSB.bin", opts->indir, opts->seg, ifo[det].name);

     if ((data = fopen(filename, "r")) != NULL) {
          // Detector position w.r.t Solar System Baricenter
          // for every datapoint
          status = fread((void *)(ifo[det].sig.DetSSB), sizeof(double), 3*sett->N, data);

          // Deterministic phase defining the position of the Earth
          // in its diurnal motion at t=0
          status = fread((void *)(&ifo[det].sig.phir), sizeof(double), 1, data);

          // Earth's axis inclination to the ecliptic at t=0
          status = fread((void *)(&ifo[det].sig.epsm), sizeof(double), 1, data);
          fclose (data);

          // printf("[%s] Using %s as ephemerids...\n", ifo[det].name, filename);
     } else {
          perror (filename);
          exit(EXIT_FAILURE);
     }

     // sincos
     ifo[det].sig.sphir = sin(ifo[det].sig.phir);
     ifo[det].sig.cphir = cos(ifo[det].sig.phir);
     ifo[det].sig.sepsm = sin(ifo[det].sig.epsm);
     ifo[det].sig.cepsm = cos(ifo[det].sig.epsm);

     sett->sepsm = ifo[det].sig.sepsm;
     sett->cepsm = ifo[det].sig.cepsm;

} // end of read detssb


/* Reads the start time of the data segment of detector det, in GPS seconds */
void read_start_time( Search_settings *sett, Command_line_opts *opts, int det )
{
     FILE *data;
     int status;
     char filename[FNAME_LENGTH];

     sprintf (filename, "%s/%03d/%s/starting_date", opts->indir, opts->seg, ifo[det].name);

     if ((data = fopen(filename, "r")) != NULL) {
          status = fscanf(data, "%lf", &ifo[det].start_time);
          fclose (data);
          printf("[%s] Starting time = %.3f\n", ifo[det].name, ifo[det].start_time);
     } else {
          perror (filename);
          exit(EXIT_FAILURE);
     }

} // end of read start time


int compared2c(const void *a, const void *b)
{
     double* da = (double*)a;
     double* db = (double*)b;

     int diff1 = (da[0] > db[0]) - (da[0] < db[0]);
     if (diff1 != 0) return diff1;
     return (da[1] > db[1]) - (da[1] < db[1]);
}


int extract_band_vlines(Search_settings *sett, Command_line_opts *opts, char *band_vl_file)
{

    int i=0, lnum, j;
    char linefile[FNAME_LENGTH], line[MAXLINE], *lfile;
    FILE *data;
    struct vlines {
        double f;
        int type;
        double offset;
        int iharm1, iharm2;
        double lwidth, rwidth;
        int det;
        int nline;
        char vfile[MAXLINE];
    } vline[MAXVFILEL];
    char line_aux[MAXVFILEL][MAXLINE];
    double fl, fr;

    glob_t globbuf;
    globbuf.gl_offs = 0;

    double fdotmin, fdotmax, Dfmaxmax;
    double normdtEmax[MAX_DETECTORS], normdEmax[MAX_DETECTORS];

    get_fdot_range(sett->fpo, sett->dt, &fdotmin, &fdotmax);

    // per-detector max line broadening due to demodulation (velocity and time derivative of DetSSB)
    for (int det=0; det < sett->nifo; det++) {

        normdEmax[det]  = 0.;
        normdtEmax[det] = 0.;

        for (i=0; i<sett->N-1; i++) {
            double dE=0, dtE=0;
            for (j=0; j<3; j++) {
                dE  += pow(fabs(ifo[det].sig.DetSSB[(i+1)*3+j] - ifo[det].sig.DetSSB[i*3+j]), 2);
                dtE += pow(fabs(ifo[det].sig.DetSSB[(i+1)*3+j]*(i+1) - ifo[det].sig.DetSSB[i*3+j]*i)*sett->dt , 2);
            }
            if (dE  > normdEmax[det])  normdEmax[det]  = dE;
            if (dtE > normdtEmax[det]) normdtEmax[det] = dtE;
        }

        printf("Detector %d: normdtEmax = %f, normdEmax = %f\n", det,
            sqrt(normdtEmax[det]), sqrt(normdEmax[det]));
    }

    // read veto files and find lines in band
    int iold = 0;
    i = 0; // veto line number (global for all veto files)

    for (int det=0; det < sett->nifo; det++) {
        // search for all veto files matching pattern <data>/lines/<det_name>lines*.csv
        sprintf(linefile, "%s/lines/%slines*.csv", opts->indir, ifo[det].name);

        printf("[%s] Looking for %s ... ", ifo[det].name, linefile);
        glob(linefile, GLOB_DOOFFS, NULL, &globbuf);
        printf("%ld files match\n", globbuf.gl_pathc);

        for (size_t ifile = 0; ifile != globbuf.gl_pathc; ++ifile){
            lfile = globbuf.gl_pathv[ifile];
            printf("   [%s] %s ", ifo[det].name, lfile);

            // Reading line data from the input file (data)
            // Columns are:
            // 1 - frequency spacing (Hz) of comb (or frequency of single line)
            // 2 - comb type (0 - singlet, 1 - comb with fixed width, 2 - comb with scaling width)
            // 3 - frequency offset of 1st visible harmonic (Hz)
            // 4 - index of first visible harmonic
            // 5 - index of last visible harmonic
            // 6 - width of left band (Hz)
            // 7 - width of right band (Hz)
            // 8 - comments

            if ((data = fopen(lfile, "r")) == NULL) {
                printf("Can't open %s \n Aborting!\n", lfile);
                exit(EXIT_FAILURE);
            }

            int nline = 0; // line number in file
            while (fgets(line, MAXLINE, data) != NULL) {
                nline++;
                // Skip comment lines beginning with '%'
                if (*line == '%') continue;
                vline[i].f       = atof(strtok(line,","));
                vline[i].type    = atoi(strtok(NULL,","));
                vline[i].offset  = atof(strtok(NULL,","));
                vline[i].iharm1  = atoi(strtok(NULL,","));
                vline[i].iharm2  = atoi(strtok(NULL,","));
                vline[i].lwidth  = atof(strtok(NULL,","));
                vline[i].rwidth  = atof(strtok(NULL,","));
                vline[i].det     = det;
                vline[i].nline = nline;
                strcpy(vline[i].vfile, lfile);
                if (++i > MAXVFILEL-1) {
                    printf("Too many lines in file %s, increase MAXVFILEL!\n", lfile);
                    exit(EXIT_FAILURE);
                }
            }

            printf("(%d data lines) \n", i-iold);
            iold = i;
            fclose(data);

        }

        globfree(&globbuf);

    } // for det


    lnum = i;

    j = 0; // index of line in band
    if (opts->narrowdown < 0.5*M_PI) j = sett->numlines_band;

    // Apply line widths
    //------------------

    for (i=0; i<lnum; i++) {

        int k;

        switch(vline[i].type) {

            // Singlet
            case 0:

            // Line width from the resampling broadening
            Dfmaxmax = 2.*fdotmax*(sett->N*sett->dt + sqrt(normdtEmax[vline[i].det])) +
                vline[i].f*sqrt(normdEmax[vline[i].det]);

            fl = vline[i].f - vline[i].lwidth - Dfmaxmax;
            fr = vline[i].f + vline[i].rwidth + Dfmaxmax;

            if (line_in_band(&fl, &fr, sett)) {
                sett->lines[j][0] = fl;
                sett->lines[j][1] = fr;
                sprintf(line_aux[j], "singlet          [l.%d]%s", vline[i].nline, vline[i].vfile);
                if (++j > MAXL-1) {printf("Too many lines, increase MAXL!\n"); exit(EXIT_FAILURE);}
            }
            break;

            // Comb with fixed width. Vetoing the band
            // [offset+index*spacing-leftwidth, offset+index*spacing+rightwidth]
            case 1:

            for (k=vline[i].iharm1; k<=vline[i].iharm2; k++) {

                double linefreq = vline[i].offset + k*vline[i].f;
                // Line width from the resampling broadening
                Dfmaxmax = 2.*fdotmax*(sett->N*sett->dt + sqrt(normdtEmax[vline[i].det])) +
                    linefreq*sqrt(normdEmax[vline[i].det]);

                fl = linefreq - vline[i].lwidth - Dfmaxmax;
                fr = linefreq + vline[i].rwidth + Dfmaxmax;

                if (line_in_band(&fl, &fr,  sett)) {
                    sett->lines[j][0] = fl;
                    sett->lines[j][1] = fr;
                    sprintf(line_aux[j], "fcomb rank %4d  [l.%d]%s", k, vline[i].nline, vline[i].vfile);
                    if (++j > MAXL-1) {printf("Too many lines, increase MAXL!\n"); exit(EXIT_FAILURE);}
                }

            }
            break;

            // Comb with scaling-width. Vetoing the band
            // [offset+index*spacing-index*leftwidth, offset+index*spacing+index*rightwidth]
            case 2:

            for (k=vline[i].iharm1; k<=vline[i].iharm2; k++) {

                double linefreq = vline[i].offset + k*vline[i].f;
                // Line width from the resampling broadening
                Dfmaxmax = 2.*fdotmax*(sett->N*sett->dt + sqrt(normdtEmax[vline[i].det])) +
                           linefreq*sqrt(normdEmax[vline[i].det]);

                fl = linefreq - k*vline[i].lwidth - Dfmaxmax;
                fr = linefreq + k*vline[i].rwidth + Dfmaxmax;

                if (line_in_band(&fl, &fr, sett)) {
                    sett->lines[j][0] = fl;
                    sett->lines[j][1] = fr;
                    sprintf(line_aux[j], "scomb rank %4d  [l.%d]%s", k, vline[i].nline, vline[i].vfile);
                    if (++j > MAXL-1) {printf("Too many lines, increase MAXL!\n"); exit(EXIT_FAILURE);}
                }
            } //k
            break;

        } // switch
    } // i

    printf("Extracted %d veto lines in band from %s [Hz, radians, line info]:\n", j-sett->numlines_band, band_vl_file);

    // write veto lines found in this band to a text file
    if ( !(data = fopen(band_vl_file, "w")) ) {
        printf("Can't open %s for writing!\n", band_vl_file);
        exit(EXIT_FAILURE);
    }
    fprintf(data, "# fl[rad]  fr[rad]   fl[Hz]  fr[Hz]   line info\n");

    sett->nvlines_all_inband = j;

    // scale veto lines to radians (narrowdown lines are already scaled)
    for (i=sett->numlines_band; i<sett->nvlines_all_inband; i++) {
        fl = sett->lines[i][0];
        fr = sett->lines[i][1];
        sett->lines[i][0] = (sett->lines[i][0] - sett->fpo)/(sett->B)*M_PI;
        sett->lines[i][1] = (sett->lines[i][1] - sett->fpo)/(sett->B)*M_PI;

        printf("   %.15f  %.15f  %.15f  %.15f  %s\n",
            sett->lines[i][0], sett->lines[i][1], fl, fr, line_aux[i]);
        fprintf(data, "   %.15f  %.15f  %.15f  %.15f  %s\n",
            sett->lines[i][0], sett->lines[i][1], fl, fr, line_aux[i]);
    }

    printf("All lines in band (incl. narrowdown, in radians):\n");
    for(i=0; i<sett->numlines_band; i++)
        printf("   %.15f %.15f (narrowdown)\n", sett->lines[i][0], sett->lines[i][1]);
    for (i=sett->numlines_band; i<sett->nvlines_all_inband; i++)
        printf("   %.15f %.15f (veto)\n", sett->lines[i][0], sett->lines[i][1]);

    // calculate veto fraction of the band (also prints it to stdout)
    double vf;
    vf = lines_veto_fraction(sett, sett->numlines_band, sett->nvlines_all_inband);
    fprintf(data, "#band_veto_fraction= %6.4f\n", vf);
    fclose(data);
    printf("Wrote veto lines in band to: %s\n", band_vl_file);

    // set number of veto lines only if veto option is given
    if (opts->veto_flag) {
         printf("Veto lines will be applied!\n");
         sett->numlines_band = sett->nvlines_all_inband;
    } else {
         printf("Veto lines WILL NOT be applied!\n");
    }

    // warnings about high veto fraction
    if( (vf > 0.9999) && opts->veto_flag) {
        printf("[Warning] This band is fully vetoed!\n");
    }
    if( (vf > 0.9) && (strcmp(opts->gtype, "allsky")==0) ){
        printf("[Warning] Band veto fraction > 90%% in allsky mode\n");
    }

    return 0;
}


int line_in_band(double* fl, double* fr, Search_settings* sett )
{
    double bs, be;        // Band start and end

    bs = sett->fpo;
    be = sett->fpo + sett->B;

    if (!(*fr < bs || *fl > be)) {
        if (*fl < bs) *fl = bs;
        if (*fr > be) *fr = be;
        return(1);
    }
    return(0);
}



void narrow_down_band(Search_settings* sett, Command_line_opts *opts)
{
    // Adding excluding ranges near the edges to the known lines list
    sett->lines[0][0] = 0;
    sett->lines[0][1] = M_PI_2 - opts->narrowdown;
    sett->lines[1][0] = M_PI_2 + opts->narrowdown;
    sett->lines[1][1] = M_PI;

    sett->numlines_band = 2;
    printf("Band is narrowed-down to [%f, %f] (narrowdown=%f)\n", sett->lines[0][1], sett->lines[1][0],
        opts->narrowdown/M_PI);

}


float lines_veto_fraction(Search_settings* sett, int lf, int le)
{

    // lf - index of first line, le - index of last line
    int i;
    double ll=0., gap=0.;
    float vf;

    // Sorting veto lines in band (1st then 2nd column)
    qsort(&sett->lines[lf], le-lf, 2*sizeof(double), compared2c);

    for (i=lf; i<le; i++) {

        // Looking for a gap between lines
        if(sett->lines[i][0] >= ll) {
            gap += sett->lines[i][0] - ll;
            ll = sett->lines[i][1];
        } else {
            if (ll < sett->lines[i][1]) ll = sett->lines[i][1];
        }
    }
    if ( ll < M_PI) gap += M_PI - ll;

    vf = (float)(M_PI-gap)/M_PI;
    printf("Band veto fraction = %6.4f\n", vf);
    return(vf);
}


void read_band_vlines(Search_settings *sett, Command_line_opts *opts, char *band_vl_file)
{
    int i, lnum;
    char line[MAXLINE];
    FILE *data;
    // band veto fraction, written by extract_band_vlines() in a comment line
    const char vf_tag[] = "#band_veto_fraction=";
    const size_t vf_taglen = sizeof(vf_tag) - 1;
    float vf = -1.;    // negative until read from the file

    i=0; // index of line in band
    if(opts->narrowdown < 0.5*M_PI) i = sett->numlines_band;

    if ((data = fopen(band_vl_file, "r")) != NULL) {
        lnum = 0;
        while (fgets(line, MAXLINE, data) != NULL) {
            lnum++;
            // Veto fraction of the band, stored after the vf_tag string
            if (!strncmp(line, vf_tag, vf_taglen)) {
                if (sscanf(line + vf_taglen, "%f", &vf) != 1) {
                    printf("Can't parse veto fraction in line %d of %s:\n%s Aborting!\n",
                        lnum, band_vl_file, line);
                    exit(EXIT_FAILURE);
                }
                continue;
            }
            // Skip other comment lines beginning with '#'
            if (*line == '#') continue;
            // Columns are whitespace separated, as written by extract_band_vlines():
            // 1 - fl [rad], 2 - fr [rad], 3 - fl [Hz], 4 - fr [Hz], 5 - line info
            if (sscanf(line, "%lf %lf", &sett->lines[i][0], &sett->lines[i][1]) != 2) {
                printf("Can't parse line %d of %s:\n%s Aborting!\n", lnum, band_vl_file, line);
                exit(EXIT_FAILURE);
            }
            if (++i > MAXL-1) {
                printf("Too many lines, increase MAXL!\n");
                exit(EXIT_FAILURE);
            }
        }
        fclose(data);
        printf("Read %d veto lines in band from %s\n", i - sett->numlines_band, band_vl_file);
    } else {
        perror (band_vl_file);
        exit(EXIT_FAILURE);
    }
    sett->nvlines_all_inband = i;

    if (vf < 0.) {
        printf("[Error] veto fraction missing in %s. Aborting!\n", band_vl_file);
        exit(EXIT_FAILURE);
    }
    printf("Band veto fraction = %6.4f\n", vf);

    printf("All lines in band (incl. narrowdown, in radians):\n");
    for(i=0; i<sett->numlines_band; i++)
        printf("   %.15f %.15f (narrowdown)\n", sett->lines[i][0], sett->lines[i][1]);
    for (i=sett->numlines_band; i<sett->nvlines_all_inband; i++)
        printf("   %.15f %.15f (veto)\n", sett->lines[i][0], sett->lines[i][1]);

    // set number of veto lines only if veto option is given
    if (opts->veto_flag) {
         printf("Veto lines will be applied!\n");
         sett->numlines_band = sett->nvlines_all_inband;
    } else {
         printf("Veto lines WILL NOT be applied!\n");
    }

    // warnings about high veto fraction
    if( (vf > 0.9999) && opts->veto_flag) {
        printf("[Warning] This band is fully vetoed!\n");
    }
    if( (vf > 0.9) && (strcmp(opts->gtype, "allsky")==0) ){
        printf("[Warning] Band veto fraction > 90%% in allsky mode\n");
    }

    return;
}
