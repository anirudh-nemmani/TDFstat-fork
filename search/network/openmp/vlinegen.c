#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <float.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <glob.h>

#include "common.h"
#include "auxi.h"
#include "struct.h"
#include "settings.h"
#include "init.h"

/*
 * vlinegen - reads a search .ini file and writes vlines.dat, the list of
 * instrumental veto lines found in the given band, plus the veto fraction.
 *
 * main() follows the same initialization sequence as the search (main.c).
 */

/* Global array of detectors (network); declared extern in struct.h,
 * defined once here since vlinegen does not link against main.c */
Detector_settings ifo[MAX_DETECTORS];


int read_lines( Search_settings *sett, Command_line_opts *opts)
{

    int i=0, lnum, j;
    char linefile[1200], line[512], *lfile;
    FILE *data;
    struct vlines {
        double f;
        int type;
        double offset;
        int iharm1, iharm2;
        double lwidth, rwidth;
        int det;
        int nline;
        char vfile[512];
    } vline[MAXVFILEL];
    char line_aux[MAXVFILEL][512];
    double fl, fr;

    glob_t globbuf;
    globbuf.gl_offs = 0;

    double fdotmin, fdotmax, Dfmaxmax;
    double normdtEmax[MAX_DETECTORS], normdEmax[MAX_DETECTORS];

    get_fdot_range(sett->fpo, sett->dt, &fdotmin, &fdotmax);

    // per-detector max line broadening due to demodulation (velocity and time derivative of DetSSB)
    for (int det=0; det < sett->nifo; det++) {

        double dEmax[3] = {0}, dtEmax[3] = {0};

        for (i=0; i<sett->N-1; i++) {
            for (j=0; j<3; j++) {
                double dE  = fabs(ifo[det].sig.DetSSB[(i+1)*3+j] - ifo[det].sig.DetSSB[i*3+j]);
                double dtE = fabs(ifo[det].sig.DetSSB[(i+1)*3+j]*(i+1) - ifo[det].sig.DetSSB[i*3+j]*i)*sett->dt;
                if (dE  > dEmax[j])  dEmax[j]  = dE;
                if (dtE > dtEmax[j]) dtEmax[j] = dtE;
            }
        }

        normdtEmax[det] = 0.;
        normdEmax[det]  = 0.;
        for (j=0; j<3; j++) {
            normdtEmax[det] += pow(dtEmax[j], 2.);
            normdEmax[det]  += pow(dEmax[j], 2.);
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
            while (fgets(line, 512, data) != NULL) {
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

    j=0; // index of line in band
    if(opts->narrowdown < 0.5*M_PI) j = sett->numlines_band;

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

    printf("%d veto lines in band [Hz, radians, line info]:\n", j-sett->numlines_band);

    // write veto lines found in this band to a text file
    sprintf(linefile, "%s/vlines_%04d.dat", opts->outdir, opts->band);
    if ( !(data = fopen(linefile, "w")) ) {
        printf("Can't open %s for writing!\n", linefile);
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
        fprintf(data, "   %.15f  %.15f  %.15f  %f  %s\n",
            sett->lines[i][0], sett->lines[i][1], fl, fr, line_aux[i]);
    }

    // calculate veto fraction of the band (also prints it to stdout)
    double vf;
    vf = lines_veto_fraction(sett, sett->numlines_band, sett->nvlines_all_inband);
    fprintf(data, "#band_veto_fraction= %6.4f\n", vf);
    fclose(data);
    printf("Wrote veto lines in band to: %s\n", linefile);

    // warnings about high veto fraction
    if( (vf > 0.999999) && opts->veto_flag) {
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


int main(int argc, char *argv[])
{

    Search_settings sett;
    Command_line_opts opts;
    struct stat buffer;
    int i;

    // Command line options (reads the .ini file given as argv[1])
    read_ini_file(&sett, &opts, argc, argv);

    // Output directory handling (create it if missing, like main.c does)
    if (stat(opts.outdir, &buffer) == -1) {
        if (errno == ENOENT) {
            if (mkdir(opts.outdir, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH) == -1) {
                perror(opts.outdir);
                return 1;
            }
        } else {
            perror(opts.outdir);
            return 1;
        }
    }

    // Detector network settings
    detectors_settings(&sett, &opts);

    // Grid data (sets sett.fftpad, needed by search_settings())
    read_grid(&sett, &opts);

    // Search settings (sets sett.N, sett.B, sett.numlines_band, ...
    // all of them used by read_lines())
    search_settings(&sett);

    // Reading the ephemerids, i.e. ifo[].sig.DetSSB needed by read_lines().
    // Only this part of init_arrays() is used here; the input time series and
    // the search work arrays are not needed to generate the veto lines.
    for (i=0; i<sett.nifo; i++)
        read_detssb(&sett, &opts, i);

    // Narrowing-down the band (excludes the edges according to
    // opts.narrowdown parameter); adds two lines
    if (opts.narrowdown < 0.5*M_PI) narrow_down_band(&sett, &opts);

    // Reading veto lines data from external files, writes vlines.dat
    printf("Reading veto files...\n");
    read_lines(&sett, &opts);

    return 0;
}
