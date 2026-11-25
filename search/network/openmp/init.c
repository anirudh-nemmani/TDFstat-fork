#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <getopt.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <time.h>
#include <float.h>

#include "common.h"
#include "auxi.h"
#include "struct.h"
#include "init.h"
#include "settings.h"
#include "../../../utils/iniparser/src/iniparser.h"

#if defined(_OPENMP)
#include <omp.h>
#endif


void read_ini_file( Search_settings *sett,
                    Command_line_opts *opts,
                    int argc,
                    char* argv[])
{

    char ini_fname[FNAME_LENGTH];
    dictionary *ini;
    int error = 0;

    if (argc > 1) {
        strcpy (ini_fname, argv[1]);
    } else if (argc > 2) {
        printf("WARNING: too many arguments, only first one is used.\n");
    } else {
        printf("ERROR: missing input file name. Call: \"<executable> search.ini\" \n");
        exit(EXIT_FAILURE);
    }

    printf ("Loading config file %s\n", ini_fname);
    if ((ini = iniparser_load(ini_fname)) == NULL) {
        perror(ini_fname);
        exit(EXIT_FAILURE);
    }
    printf("-- INI file contents --\n");
    iniparser_dump(ini, stdout);
    printf("-----------------------\n");

    // directory containing input data
    opts->indir = iniparser_getstring(ini, "search:indir", NULL);
    // output directory
    opts->outdir = iniparser_getstring(ini, "search:outdir", ".");
    // band number
    opts->band = iniparser_getint(ini, "search:band", 0);
    // time segment number
    opts->seg = iniparser_getint(ini, "search:seg", 0);
    // hemisphere 1,2 or 0 for both
    opts->hemi = iniparser_getint(ini, "search:hemi", 0);
    // candidate fstatistic threshold
    opts->thr = iniparser_getdouble(ini, "search:thr", 20.);
    // number of days (integer)
    sett->nod = iniparser_getint(ini, "search:nod", 0);
    // sampling interval of the input time series [seconds] (double)
    sett->dt = iniparser_getdouble(ini, "search:dt", -1.);
    // bands overlap [0-1.] defines band base frequency, fpo
    opts->overlap = iniparser_getdouble(ini, "search:overlap", 0.);
    // [0-0.5] range of frequencies to output; 0 is nothing, 0.5 is full band,
    // -1 means automatic for given overlap
    opts->narrowdown = iniparser_getdouble(ini, "search:narrowdown", -1.);
    // use data from subset of detectors only (default is to use all available)
    opts->usedet = iniparser_getstring(ini, "search:usedet", "");
    // path to the grid file (absolute or relative)
    opts->grid_file = iniparser_getstring(ini, "search:grid_file", "");
    // name of the range file to read
    // opts->range_file = iniparser_getstring(ini, "search:range_file", "");
    // name of the range file to dump and exit
    // opts->dump_range_file = iniparser_getstring(ini, "search:dump_range_file", "");
    // name of the file with signall to add
    // TODO: single line signals to enable multiple signals, more flexible gsize
    opts->addsig = iniparser_getstring(ini, "search:addsig", "");
    // optional label of output files
    opts->label = iniparser_getstring(ini, "search:label", "");
    if (strlen(opts->label)) {
        const char *tmp = opts->label;
        opts->label = malloc(strlen(tmp) + 2);
        sprintf((char*)opts->label, "_%s", tmp);
    }
    // runtime modifiers, supported values are: {read_O3}
    opts->mods = iniparser_getstring(ini, "search:mods", "");

    // fstatistics normalization; if NULL white noise is assumed
    // currently only NULL and old blocks_avg are implemented
    opts->fstat_norm = iniparser_getstring(ini, "search:fstat_norm", NULL);

    //   FLAGS
    // [0, 1] (in future other values can select different subsets of lines
    opts->veto_flag = iniparser_getint(ini, "search:veto_flag", 0);
    // [0, 1] generate checkpoint files at each new sky position
    opts->checkp_flag = iniparser_getint(ini, "search:checkp_flag", 0);

    // grid specification
    opts->gtype = iniparser_getstring(ini, "grid:gtype", "allsky");
    opts->gcenter = iniparser_getstring(ini, "grid:gcenter", NULL);
    opts->gsizes = iniparser_getstring(ini, "grid:gsizes", NULL);
    opts->gsteps = iniparser_getstring(ini, "grid:gsteps", NULL);

    // injections configuration
    opts->use_amp = iniparser_getint(ini, "injections:use_amp", 1);
    opts->inj_mode = iniparser_getstring(ini, "injections:inj_mode", NULL);
    opts->inj_file = iniparser_getstring(ini, "injections:inj_file", NULL);
    opts->inj_dict = iniparser_getstring(ini, "injections:inj_dict", NULL);
    // empty values (key = ) are NULL, like missing keys
    if (opts->inj_mode && !opts->inj_mode[0]) opts->inj_mode = NULL;
    if (opts->inj_file && !opts->inj_file[0]) opts->inj_file = NULL;
    if (opts->inj_dict && !opts->inj_dict[0]) opts->inj_dict = NULL;

    // injection options
    if (opts->inj_mode) {
        if (strcmp(opts->inj_mode, "serial") != 0 && strcmp(opts->inj_mode, "multi") != 0) {
            printf("[ERROR] invalid inj_mode: %s. Must be 'serial' or 'multi'.\n", opts->inj_mode);
            exit(EXIT_FAILURE);
        }
        if (!opts->inj_file && !opts->inj_dict) {
            printf("[ERROR] inj_mode = %s needs inj_file or inj_dict\n", opts->inj_mode);
            exit(EXIT_FAILURE);
        }
        if (opts->inj_file && opts->inj_dict) {
            printf("[ERROR] inj_file and inj_dict cannot be used together. Please specify only one.\n");
            exit(EXIT_FAILURE);
        }
        if (opts->inj_file) {
            // check if the file exists
            FILE *file = fopen(opts->inj_file, "r");
            if (file) {
                fclose(file);
            } else {
                printf("[ERROR] inj_file can't be opened: %s\n", opts->inj_file);
                exit(EXIT_FAILURE);
            }
        }
    }

    // gcenter = injection (not for allsky, which ignores gcenter): a grid center for every injection,
    // so only with inj_mode = serial
    if (strcmp(opts->gtype, "allsky") != 0 && opts->gcenter && strcmp(opts->gcenter, "injection") == 0 &&
        !(opts->inj_mode && strcmp(opts->inj_mode, "serial") == 0)) {
        printf("[ERROR] gcenter = injection needs inj_mode = serial\n");
        exit(EXIT_FAILURE);
    }

    // various checks
    if (! opts->indir) {
        error = 1; printf("[ERROR] missing indir !\n");
    }
    if (opts->band < 1) {
        error = 1; printf("[ERROR] missing band !\n");
    }
    if (opts->seg < 1) {
        error = 1; printf("[ERROR] missing segment !\n");
    }
    if (sett->nod < 1) {
        error = 1; printf("[ERROR] missing nod !\n");
    }
    if (sett->dt < 0.) {
        error = 1; printf("[ERROR] missing dt !\n");
    }

    if (error == 1) exit(EXIT_FAILURE);

    // test overlap and narrowdown
    // todo: NAV
    if (opts->narrowdown < 0.) {
        opts->narrowdown = (1. - opts->overlap)/2. ;
    }
    opts->narrowdown *= M_PI;

    //sett->fpo = 10. + (1. - opts->overlap)*opts->band*(0.5/sett->dt);
    sett->fpo = fpo(opts->band, opts->overlap, sett->dt);
} // read_ini_file


/* Generate grid from the M matrix (grid.bin) */

void read_grid( Search_settings *sett, Command_line_opts *opts )
{

     //sett->M = (double *) calloc (16, sizeof (double));

     FILE *data;
     char filename[FNAME_LENGTH];
     int i;

     // In case when usedet option is used for one detector
     // i.e. opts->usedet has a length of 2 (e.g. H1 or V1),
     // or data dir contains only one detector,
     // read grid.bin from this detector subdirectory
     // (see detectors_settings() in settings.c for details)

     if ( opts->grid_file[0] == '/' ) {
          // absolute path specified
          sprintf (filename, "%s", opts->grid_file);
     } else if ( strlen(opts->grid_file) > 0 ) {
          // relative path specified
          sprintf (filename, "%s/%s", opts->indir, opts->grid_file);
     } else if (strlen(opts->usedet)==2 && sett->nifo==1) {
          // default path for one detector
          sprintf (filename, "%s/%03d/%s/grid.bin", opts->indir, opts->seg, ifo[0].name);
     } else {
          // default path for network of detectors
          char dnet_str[MAX_DETECTORS*DETNAME_LENGTH]="";
          for(i=0; i<sett->nifo; i++)
               strcat(dnet_str, ifo[i].name);
          sprintf (filename, "%s/%03d/grids/grid_%03d_%04d_%s.bin", opts->indir, opts->seg, opts->seg, opts->band, dnet_str);
     }

     if ((data=fopen (filename, "r")) != NULL) {
          printf("grid_file = %s\n", filename);
          fread ((void *)&sett->fftpad, sizeof (int), 1, data);
          printf("fftpad from the grid file: %d\n", sett->fftpad);
          // M: vector of 16 components consisting of 4 rows
          // of 4x4 grid-generating matrix
          fread ((void *)sett->M, sizeof (double), 16, data);
          printf("Grid generating matrix M=\n");
          for (i=0; i<4; i++)
               printf("  %12.6e  %12.6e  %12.6e  %12.6e\n",
                    sett->M[i], sett->M[4+i], sett->M[8+i], sett->M[12+i]);
          fclose (data);
     } else {
          perror (filename);
          exit(EXIT_FAILURE);
     }

} // end of read grid


/* Array initialization */

/* Allocates the per-detector work arrays of the search (amplitude modulation
 * and time shift arrays). Not needed by tools that only read the input data.
 */
void alloc_sig_array( Search_settings *sett, int det )
{

     ifo[det].sig.xDatma = fftw_malloc(sett->N*sizeof(complex double));
     ifo[det].sig.xDatmb = fftw_malloc(sett->N*sizeof(complex double));

     ifo[det].sig.aa = (double *) calloc(sett->N, sizeof(double));
     ifo[det].sig.bb = (double *) calloc(sett->N, sizeof(double));

     ifo[det].sig.shft = (double *) calloc(sett->N, sizeof(double));
     ifo[det].sig.shftf = (double *) calloc(sett->N, sizeof(double));

} // end of alloc sig array


/* Auxiliary arrays, Earth's rotation */

void init_aux_array( Search_settings *sett, Aux_arrays *aux_arr )
{

     int i;
     double omrt;

     aux_arr->t2 = (double *) calloc(sett->N, sizeof (double));
     aux_arr->cosmodf = (double *) calloc(sett->N, sizeof (double));
     aux_arr->sinmodf = (double *) calloc(sett->N, sizeof (double));

     for (i=0; i<sett->N; i++) {
          omrt = (sett->omr)*i;     // Earth angular velocity * dt * i
          aux_arr->t2[i] = sqr((double)i);
          aux_arr->cosmodf[i] = cos(omrt);
          aux_arr->sinmodf[i] = sin(omrt);
     }

} // end of init aux array


/* Allocates and initializes the data, detector ephemeris and F-statistic
 * arrays, by calling the individual readers (common.c) and allocators for
 * every detector of the network.
 */

void init_arrays( Search_settings *sett,
                  Command_line_opts *opts,
                  Aux_arrays *aux_arr )
{

     int i;

     for (i=0; i<sett->nifo; i++) {

          // Input time-domain data handling
          read_xdat(sett, opts, i);

          // Ephemeris file handling
          read_detssb(sett, opts, i);

          // Start time reading
          read_start_time(sett, opts, i);

          alloc_sig_array(sett, i);

     } // end loop for detectors

     // Safe check for the start time
     double st_temp = ifo[0].start_time;
     for (i=1; i<sett->nifo; i++) {
          if (ifo[i].start_time != st_temp) {
              printf("Start time doesn't match between detectors %s and %s. Aborting...\n", ifo[0].name, ifo[i].name);
              exit(EXIT_FAILURE);
          }
     }


     // Check if the ephemerids have the same epsm parameter
     for(i=1; i<sett->nifo; i++) {
          if(!(ifo[i-1].sig.sepsm == ifo[i].sig.sepsm)) {
               printf("The parameter epsm (DetSSB.bin) differs for detectors %s and %s. Aborting...\n",
                    ifo[i-1].name, ifo[i].name);
               exit(EXIT_FAILURE);
          }
     }

     // if all is well with epsm, take the first value
     sett->sepsm = ifo[0].sig.sepsm;
     sett->cepsm = ifo[0].sig.cepsm;

     init_aux_array(sett, aux_arr);

} // end of init arrays


/* Allocates and reads the injection file or string in INI file to
 * be later used for add_injection function. Blank lines and lines
 * starting with '#' are ignored and returns the number of injections.
 * Format: id t_ref amp f fdot ra dec iota psi phase
 * WIP - When we use TOML format for INI files, we need to adjust this
 * for the TOML table format.
 */

int read_injection( Command_line_opts *opts,
                    Injection_params **inj)
{
    FILE *data = NULL;
    char line[MAXLINE];
    char c;
    int ninj = 0;
    Injection_params ip;

    *inj = NULL;

    // open both sources as a FILE*
    if (opts->inj_file)
        data = fopen(opts->inj_file, "r");
    else if (opts->inj_dict)
        // WIP: When TOML INI parser comes live, we have to change fmemopen with
        // something that fills Injection_params from TOML table.
        data = fmemopen((void *)opts->inj_dict, strlen(opts->inj_dict), "r");

    while (data && fgets(line, sizeof(line), data)) {
        // skip blank lines and comments: first non-whitespace char missing or '#'
        if (sscanf(line, " %c", &c) != 1 || c == '#') continue;

        if (sscanf(line, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                &ip.id, &ip.t_ref, &ip.amp, &ip.f, &ip.fdot,
                &ip.ra, &ip.dec, &ip.iota, &ip.psi, &ip.phase) != 10) {
            printf("[ERROR] read_injection(): malformed injection line:\n   '%s'\n", line);
            exit(EXIT_FAILURE);
        }
        *inj = realloc(*inj, (ninj + 1) * sizeof(Injection_params));
        (*inj)[ninj++] = ip;
    }
    if (data) fclose(data);

    if (ninj == 0) {
        printf("[ERROR] read_injection(): inj_mode is set but no injections were found in inj_file or inj_dict.\n"
               "        Set inj_mode empty to run without injections.\n");
        exit(EXIT_FAILURE);
    }

    printf("read_injection(): %d injection(s) read\n", ninj);
    for (int i = 0; i < ninj; i++)
        printf("   id=%d t_ref=%f %s=%le f=%le fdot=%le ra=%le dec=%le iota=%le psi=%le phase=%le\n",
            (*inj)[i].id, (*inj)[i].t_ref, opts->use_amp ? "amp" : "snr", (*inj)[i].amp, (*inj)[i].f, (*inj)[i].fdot,
            (*inj)[i].ra, (*inj)[i].dec, (*inj)[i].iota, (*inj)[i].psi, (*inj)[i].phase);

    return ninj;
} // end of read_injection


/* Generates the injection and adds them to the xDat array in the
 * detector structures. If the injection mode is "multi", it adds
 * all the injections to the xDat array. If it's "serial" mode, it
 * adds one by one based on the index idx from the main.c loop.
 */

void add_injection( Search_settings *sett,
                    Command_line_opts *opts,
                    Aux_arrays *aux_arr,
                    Injection_params *inj,
                    int ninj,
                    int inj_idx)
{
    int i, j, n, k, nadd = 0;
    double sum, h0, delta;
    double sinalpha, cosalpha, sindelta, cosdelta, phaseadd, shiftadd;
    double phi, psi, cosi, cosip, amplit[4];
    double nSource[3], freqo[2];

    // Storing the original noise data
    // Allocated only when injection mode is "serial"
    static float *xDat_orig[MAX_DETECTORS];

    if (!strcmp(opts->inj_mode, "serial")) {
        for (n=0; n<sett->nifo; n++) {
            if (!xDat_orig[n]) {
                // first call: store the original data
                xDat_orig[n] = (float *)malloc(sett->N*sizeof(float));
                memcpy(xDat_orig[n], ifo[n].sig.xDat, sett->N*sizeof(float));
            } else {
                // next calls: start again from the original data
                memcpy(ifo[n].sig.xDat, xDat_orig[n], sett->N*sizeof(float));
            }
        }
    }

    // Allocate arrays for added injection, for each detector
    double **injadd = malloc((sett->nifo)*sizeof(double *));
    for (n=0; n<sett->nifo; n++)
        injadd[n] = malloc((sett->N)*sizeof(double));

    // Allocation for saving the injected signal if opts->mods has "save_inj"
    int save_inj = (opts->mods && strstr(opts->mods, "save_inj") != NULL);
    double **injtot = NULL;
    if (save_inj && !strcmp(opts->inj_mode, "multi")) {
        injtot = malloc((sett->nifo)*sizeof(double *));
        for (n=0; n<sett->nifo; n++)
            injtot[n] = calloc(sett->N, sizeof(double));
    }

    // Add all injections and skip if the injection mode is "serial"
    for (k=0; k<ninj; k++) {

        if (!strcmp(opts->inj_mode, "serial") && k != inj_idx) continue;

        Injection_params *ip = &inj[k];
        sum = 0.;

        // Offset (in samples) between the reference time and the start of the current segment
        delta = (ifo[0].start_time - ip->t_ref)/sett->dt;

        // Converting the frequency into dimensionless units
        freqo[0] = 2 * M_PI * ip->f * sett->dt;
        freqo[1] = M_PI * ip->fdot * sett->dt * sett->dt;

        // Shift the frequency to the current segment based on spindown from the reference time
        freqo[0] += 2.*freqo[1]*delta;

        // Check if the signal is in band
        if( freqo[0]/(2*M_PI*sett->dt) < sett->fpo || freqo[0]/(2*M_PI*sett->dt) > sett->fpo + sett->B ) {
            printf("[WARNING] injection id=%d out of band f=%le s=%le\n",
                   ip->id, freqo[0]/(2*M_PI*sett->dt), ip->fdot);
            continue;
        }

        // Calculation of sin alpha, cos alpha, sin delta, cos delta of the signal.
        // Check Eq. 18 of Phys. Rev. D 58, 063001 1998
        sinalpha = sin(ip->ra);
        cosalpha = cos(ip->ra);
        sindelta = sin(ip->dec);
        cosdelta = cos(ip->dec);

        // Calculation of four amplitudes from polarization, phase and inclination
        // Check Eq. 32 - 35 of Phys. Rev. D 58, 063001 1998
        cosi = cos(ip->iota);
        psi = ip->psi;
        phi = ip->phase;
        cosip = (1. + cosi*cosi)/2.;

        amplit[0] = cos(2.*psi)*cosip*cos(phi) - sin(2.*psi)*cosi*sin(phi);
        amplit[1] = sin(2.*psi)*cosip*cos(phi) + cos(2.*psi)*cosi*sin(phi);
        amplit[2] = -cos(2.*psi)*cosip*sin(phi) - sin(2.*psi)*cosi*cos(phi);
        amplit[3] = -sin(2.*psi)*cosip*sin(phi) + cos(2.*psi)*cosi*cos(phi);

        // To keep coherent phase between time segments
        double phaseshift = freqo[0]*delta - freqo[1]*pow(delta, 2);

        // Loop for each detector - sum calculations
        for (n=0; n<sett->nifo; n++) {

            modvir(sinalpha, cosalpha, sindelta, cosdelta, sett->N, &ifo[n], aux_arr);

            nSource[0] = cosalpha*cosdelta;
            nSource[1] = sinalpha*cosdelta;
            nSource[2] = sindelta;

            for (i=0; i<sett->N; i++) {

                // Gaps in the data (xDat = 0): no signal, not counted in the SNR
                if (fabs(ifo[n].sig.xDat[i]) <= FLT_MIN) {
                    injadd[n][i] = 0.;
                    continue;
                }

                // Calculation of n*r / c*delta_t vector
                shiftadd = 0.;
                for (j=0; j<3; j++)
                    shiftadd += nSource[j]*ifo[n].sig.DetSSB[i*3+j];

                // Phase (no-approximations)
                // Matches with the solution in LAL Pulsar
                // Refer Phys. Rev. D 59, 063003 1999
                phaseadd = freqo[0]*(i + shiftadd) + freqo[1]*(2.*shiftadd*i + aux_arr->t2[i] + shiftadd*shiftadd) + phaseshift;
                // We now dephase this to shift the signal within the bandwidth
                phaseadd -= 2*M_PI*sett->dt*sett->fpo*i;

                // The whole signal with 4 amplitudes and modulations
                injadd[n][i] = amplit[0]*(ifo[n].sig.aa[i])*cos(phaseadd)
                    + amplit[1]*(ifo[n].sig.bb[i])*cos(phaseadd)
                    + amplit[2]*(ifo[n].sig.aa[i])*sin(phaseadd)
                    + amplit[3]*(ifo[n].sig.bb[i])*sin(phaseadd);

                // Noise-weighted sum over signal (for the network SNR)
                sum += pow(injadd[n][i], 2.)/ifo[n].sig.sig2;

            } // data loop
        } // detector loop

        // Signal amplitude h0: 3rd column is amplitude (use_amp=1) or network SNR (use_amp=0)
        // network SNR^2 = h0^2 * sum_n sum_i s_n[i]^2 / sig2_n  (Gaussian noise with fixed sigma per detector)
        if (opts->use_amp) {
            h0 = ip->amp;
        } else {
            h0 = ip->amp / sqrt(sum);
        }
        printf("add_injection(): injection id=%d h0=%le\n", ip->id, h0);

        // Loop for each detector - scaling by h0 and adding signal to data (point by point)
        for (n=0; n<sett->nifo; n++) {
            for (i=0; i<sett->N; i++) {
                injadd[n][i] *= h0;
                // Adding the signal to the data vector
                ifo[n].sig.xDat[i] += injadd[n][i];
                if (injtot)
                    injtot[n][i] += injadd[n][i];
            } // data loop
        } // detector loop
        nadd++;

    } // injections loop

    // Saving the injected signal, one file per detector (skipped if all injections were out of band)
    // suffix as the output dataset name: inj_<id> (serial) or inj_0 (multi)
    if (save_inj && nadd > 0) {
        char filename[FNAME_LENGTH];
        int inj_id = strcmp(opts->inj_mode, "serial") ? 0 : inj[inj_idx].id;

        snprintf(filename, FNAME_LENGTH, "%s/injections", opts->outdir);
        mkdir(filename, 0755);

        for (n=0; n<sett->nifo; n++) {
            snprintf(filename, FNAME_LENGTH, "%s/injections/signal_%03d_%04d_%s_inj_%d%s",
                opts->outdir, opts->seg, opts->band, ifo[n].name, inj_id, opts->label);
            FILE *out = fopen(filename, "w");
            if (out == NULL) {
                perror(filename);
                exit(EXIT_FAILURE);
            }

            for (i=0; i<sett->N; i++) {
                fprintf(out, "%e\n", injtot ? injtot[n][i] : injadd[n][i]);
            }

            fclose(out);
        } // detector loop
    } // writing signal

    if (injtot) {
        for (n=0; n<sett->nifo; n++)
            free(injtot[n]);
        free(injtot);
    }

    // Free auxiliary 2d array
    for (n=0; n<sett->nifo; n++)
        free(injadd[n]);
    free(injadd);

} // add_injection()



/* Search range */

void set_search_range( Search_settings *sett,
		               Command_line_opts *opts,
                       Search_range *s_range,
                       Injection_params *inj)
{
    /*
    Sets the grid range in code units.
    Input: opts, sett, inj (current injection for gcenter = injection, otherwise NULL)
    Output: s_range

    The grid is formed by points in 4D search parameters space:
            frequency, spindown, two sky positions (m, n) which map to (ra, dec).
    The reference grid is standard, allsky, "integer" grid with integer node coordinates and spacing of 1.
    For allsky_range and directed_range fractional grid points are allowed.
    Grid range is specified by three 4D vectors:
        - gcenter: middle point of the grid in physical units (f [Hz], fdot [Hz/s], ra and dec [radians])
        - gsizes: grid range from gcenter, so the grid span is [gcenter-gsizes, gcenter+gsizes];
                  gsizes units = grid spacing of the reference, allsky grid (1)
        - gsteps: spacing of grid points in units of the reference grid spacing (1)
    Example for 1D: gsize=2.5 and gstep=0.5 results in range [gcenter-5, gcenter+5] with spacing 0.5,
                which means 19 grid points with density twice the density of allsky grid.

    There are 3 grid types (gtype):
        * allsky: the whole sky is searched; coordinates of grid points are integer with spacing of 1;
                  gcenter, gsizes, gsteps are ignored; gsteps=1.
        * allsky_range: gcenter is shifted to the nearest integer grid point;
                        fractional grid range is +/- gsizes from gcenter;
                        fractional grid spacing is set to gsteps.
        * directed_range: gcenter is exactly at the specified coordinates;
                          gsizes and gsteps can be frectional like for allsky_range.
    */

    float fr,
          sr, sr_gsize, sr_step,
          alphar, deltar, nr_gsize, nr_step, mr_gsize, mr_step;
    int fr_gsize, fr_step, fmid;
    float smid, nmid, mmid;
    double sgnlol[4];

    // Hemispheres (with respect to the ecliptic)
    if(opts->hemi) {
        s_range->pmr[0] = opts->hemi;
        s_range->pmr[1] = opts->hemi;
    } else {
        s_range->pmr[0] = 1;
        s_range->pmr[1] = 2;
    }
    // default grid steps
    s_range->mstep = 1.;
    s_range->nstep = 1.;
    s_range->sstep = 1.;

    if ( strcmp(opts->gtype, "allsky")==0 ) {
        // Maximum grid ranges given grid matrix M
        gridr(sett->M, s_range->spndr, s_range->nr, s_range->mr, sett->oms, sett->Smax);
        s_range->fr[0] = sett->nmin;
        s_range->fr[1] = sett->nmax;

        printf("[set_search_range] gtype = allsky\n");
        printf("[set_search_range] Hemispheres = %d %d\n", s_range->pmr[0], s_range->pmr[1]);
        printf("[set_search_range] Maximum grid ranges:\n");
        printf("---------------------------------------\n");
        printf("   |          min              max\n");
        printf("---------------------------------------\n");
        printf("mm | %12g     %12g\n", s_range->mr[0], s_range->mr[1]);
        printf("nn | %12g     %12g\n", s_range->nr[0], s_range->nr[1]);
        printf("s  | %12g     %12g\n", s_range->spndr[0], s_range->spndr[1]);
        printf("f  | %12d     %12d\n", s_range->fr[0], s_range->fr[1]);
        printf("---------------------------------------\n\n");

    } else {
        // setup grid center
        if (strcmp(opts->gcenter, "injection")==0) {
            if (inj == NULL) {
                printf("[ERROR] gcenter = injection but no injection given\n");
                exit(EXIT_FAILURE);
            }
            // frequency at the start of the segment (spindown from t_ref, as in add_injection)
            fr = inj->f + inj->fdot*(ifo[0].start_time - inj->t_ref);
            sr = inj->fdot;
            alphar = inj->ra;
            deltar = inj->dec;
            printf("Setting grid center at injection id=%d\n", inj->id);
            if (fr < sett->fpo || fr > sett->fpo + sett->B) {
                printf("[ERROR] injection id=%d out of band f=%g, grid center outside the band\n", inj->id, fr);
                exit(EXIT_FAILURE);
            }
        } else {
            printf("Reading grid center\n");
            if (sscanf(opts->gcenter, "%f %f %f %f", &fr, &sr, &alphar, &deltar) != 4) {
                printf("[ERROR] gcenter = %s: must be 'f fdot ra dec' or 'injection'\n", opts->gcenter);
                exit(EXIT_FAILURE);
            }
        }

        // setup grid sizes
        if (opts->gsizes) {
            sscanf(opts->gsizes, "%d %f %f %f", &fr_gsize, &sr_gsize, &mr_gsize, &nr_gsize);
        } else {
            printf("[ERROR] missing gsizes!\n");
            exit(EXIT_FAILURE);
        }

        // setup grid steps
        if (opts->gsteps) {
            sscanf(opts->gsteps, "%d %f %f %f", &fr_step, &sr_step, &mr_step, &nr_step);
        } else {
            printf("[ERROR] missing gsteps!\n");
            exit(EXIT_FAILURE);
        }

        fr_step = 1; // fr_step is always 1 (just FFT bin numbers)

        // convert fr, sr, alphar, deltar to fractional code units

        // first convert from physical to dimensionless units
        //sgnlol[0] = (fr - sett->fpo)/sett->B * M_PI;
        sgnlol[0] = 0.;
        sgnlol[1] = M_PI * sr * sett->dt * sett->dt;

        double cof = sett->oms; //+ sgnlol[0];
        double be[2];
        s_range->pmr[0] = ast2lin((double)alphar, (double)deltar, C_EPSMA, be);
        s_range->pmr[1] = s_range->pmr[0];
        sgnlol[2] = be[0]*cof;
        sgnlol[3] = be[1]*cof;

        // dimensionless to code units conversion
        float *sgnl_code = (float *) calloc(4, sizeof(float));
        dimless_to_code(sett->M, sgnlol, sgnl_code);

        // reverse conversions from jobcore()
        FLOAT_TYPE sgnl0;
        double het1, het0;
        het0 = fmod(sgnl_code[2]*sett->M[8] + sgnl_code[3]*sett->M[12], sett->M[0]);
        het1 = fmod(sgnl_code[1]*sett->M[4], sett->M[0]);
        if(het1<0) het1 += sett->M[0];
        sgnl0 = het0 + het1;

        sgnl_code[0] = (fr - sett->fpo)/sett->B * sett->nfftf/2 - sgnl0*sett->nfftf/(2*M_PI);
        //printf("di(sgnl0)=%g\n", sgnl0*sett->nfftf/(2*M_PI));

        // Set the steps
        s_range->mstep = mr_step;
        s_range->nstep = nr_step;
        s_range->sstep = sr_step;

        // round gsizes to integer multiples of steps
        fr_gsize = round(fr_gsize);
        sr_gsize = round(sr_gsize/sr_step)*sr_step;
        nr_gsize = round(nr_gsize/nr_step)*nr_step;
        mr_gsize = round(mr_gsize/mr_step)*mr_step;

        if ( strcmp(opts->gtype, "allsky_range") == 0) {

            // nearest "integer" grid point to alphar, deltar, sr, fr.
            fmid = round(sgnl_code[0]);
            smid = round(sgnl_code[1]);
            nmid = round(sgnl_code[2]);
            mmid = round(sgnl_code[3]);

        } else if (strcmp(opts->gtype, "directed_range") == 0) {

            // exact position: alphar, deltar, sr, fr.
            fmid = sgnl_code[0];
            smid = sgnl_code[1];
            nmid = sgnl_code[2];
            mmid = sgnl_code[3];

        } else {
            printf("Unknown grid type %s \n", opts->gtype);
            exit(EXIT_FAILURE);
        }

        s_range->fr[0] = fmid - fr_gsize;
        s_range->fr[1] = fmid + fr_gsize;
        s_range->spndr[0] = smid - sr_gsize;
        s_range->spndr[1] = smid + sr_gsize;
        s_range->nr[0] = nmid - nr_gsize;
        s_range->nr[1] = nmid + nr_gsize;
        s_range->mr[0] = mmid - mr_gsize;
        s_range->mr[1] = mmid + mr_gsize;

        free(sgnl_code);

        // TODO: check for out of scale ranges

        printf("[set_search_range] gtype = %s\n", opts->gtype);
        printf("[set_search_range] gcenter: f=%g [Hz]  fdot=%g [Hz/s]  ra=%g [rad]  de=%g [rad]\n", fr, sr, alphar, deltar);
        printf("[set_search_range] hemisphere = %d\n", s_range->pmr[0]);
        printf("[set_search_range] Actual, computed grid ranges:\n");
        printf("-----------------------------------------------------------------------\n");
        printf("   |            min           center              max             step \n");
        printf("-----------------------------------------------------------------------\n");
        printf("mm |   %12g     %12g     %12g     %12g\n", s_range->mr[0], mmid, s_range->mr[1], s_range->mstep);
        printf("nn |   %12g     %12g     %12g     %12g\n", s_range->nr[0], nmid, s_range->nr[1], s_range->nstep);
        printf("s  |   %12g     %12g     %12g     %12g\n", s_range->spndr[0], smid, s_range->spndr[1], s_range->sstep);
        printf("f  |   %12d     %12d     %12d     %12d\n", s_range->fr[0], fmid, s_range->fr[1], 1);
        printf("-----------------------------------------------------------------------\n\n");

    }
} // end of set_search_range


/* FFT Plans	 */

void plan_fftw( Search_settings *sett,
                Command_line_opts *opts,
                FFTW_plans *plans,
                FFTW_arrays *fftw_arr,
                Aux_arrays *aux_arr)
{

#if defined(_OPENMP)
     fftw_init_threads();
     fftw_plan_with_nthreads(omp_get_max_threads());
#ifdef COMP_FLOAT
     fftwf_init_threads();
     FFTW_PRE(_plan_with_nthreads)(omp_get_max_threads());
#endif
#endif

//#define FFTW_RIGOR FFTW_ESTIMATE
#define FFTW_RIGOR FFTW_MEASURE

     // arrays xa,xb are used for in-place interpolation,
     // thus their length is max{fftpad*nfft, Ninterp}
     fftw_arr->arr_len = (sett->nfftf > sett->Ninterp ? sett->nfftf : sett->Ninterp);

     fftw_arr->xa = (fftw_complex *)fftw_malloc(fftw_arr->arr_len*sizeof(fftw_complex));
     fftw_arr->xb = (fftw_complex *)fftw_malloc(fftw_arr->arr_len*sizeof(fftw_complex));

     printf("[fft plans] arr_len=%d   nfft=%d   fftpad=%d  Ninterp=%d\n",
            fftw_arr->arr_len, sett->nfft, sett->fftpad, sett->Ninterp);

     printf("double plans...");
     // Change FFTW_MEASURE to FFTW_PATIENT for more optimized plan
     plans->pl_int = fftw_plan_dft_1d(sett->nfft, fftw_arr->xa, fftw_arr->xa,
                                      FFTW_FORWARD, FFTW_RIGOR);
     plans->pl_inv = fftw_plan_dft_1d(sett->Ninterp, fftw_arr->xa, fftw_arr->xa,
                                      FFTW_BACKWARD, FFTW_RIGOR);
     printf("done\n");

     printf("float plans...");
     fftw_arr->fxa = (FFTW_PRE(_complex) *)FFTW_PRE(_malloc)(sett->nfftf*sizeof(FFTW_PRE(_complex)));
     fftw_arr->fxb = (FFTW_PRE(_complex) *)FFTW_PRE(_malloc)(sett->nfftf*sizeof(FFTW_PRE(_complex)));

     plans->plan = FFTW_PRE(_plan_dft_1d)(sett->nfftf, fftw_arr->fxa, fftw_arr->fxa,
				        FFTW_FORWARD, FFTW_RIGOR);
     printf("done\n");

} // end of FFT plans


/* Checkpointing
it is kept here but we don't use it anymore */

void read_checkpoints( Command_line_opts *opts,
                       Search_range *s_range,
                       int *FNum)
{

     if(opts->checkp_flag) {

          // filename of checkpoint state file, depending on the hemisphere
          if(opts->hemi){
               sprintf(opts->state_file, "state_%03d_%04d%s_%d.dat",
                       opts->seg, opts->band, opts->label, opts->hemi);
          } else {
               sprintf(opts->state_file, "state_%03d_%04d%s.dat",
                       opts->seg, opts->band, opts->label);
          }

          FILE *state;
          if((state = fopen(opts->state_file, "r")) != NULL) {

               // Scan the state file to get last recorded parameters
               if((fscanf(state, "%d %f %f %f %d", &s_range->pst, &s_range->mst,
                                                   &s_range->nst, &s_range->sst, FNum)) == EOF) {

                    // This means that state file is empty (=end of the calculations)
                    fprintf (stderr, "State file empty: nothing to do...\n");
                    fclose (state);
                    exit(EXIT_FAILURE);
               }

               fclose (state);

               // No state file - start from the beginning
          } else {
               s_range->pst = s_range->pmr[0];
               s_range->mst = s_range->mr[0];
               s_range->nst = s_range->nr[0];
               s_range->sst = s_range->spndr[0];
               *FNum = 0;
          } // if state

     } else {  // no checkpointing
          s_range->pst = s_range->pmr[0];
          s_range->mst = s_range->mr[0];
          s_range->nst = s_range->nr[0];
          s_range->sst = s_range->spndr[0];
          *FNum = 0;
     } // if checkp_flag

} // end reading checkpoints


  /* Cleanup & memory free  */

void cleanup( Search_settings *sett,
              Command_line_opts *opts,
              Search_range *s_range,
              FFTW_plans *plans,
              FFTW_arrays *fftw_arr,
              Aux_arrays *aux)
{

     int i;

     for(i=0; i<sett->nifo; i++) {
          free(ifo[i].sig.xDat);
          free(ifo[i].sig.xDatma);
          free(ifo[i].sig.xDatmb);
          free(ifo[i].sig.DetSSB);
          free(ifo[i].sig.aa);
          free(ifo[i].sig.bb);
          free(ifo[i].sig.shftf);
          free(ifo[i].sig.shft);
     }

     free(aux->sinmodf);
     free(aux->cosmodf);
     free(aux->t2);

     fftw_free(fftw_arr->xa);
     fftw_free(fftw_arr->xb);

     //free(sett->M);

     FFTW_PRE(_destroy_plan)(plans->plan);
     fftw_destroy_plan(plans->pl_int);
     fftw_destroy_plan(plans->pl_inv);
     /*
     fftw_destroy_plan(plans->pl_int2);
     fftw_destroy_plan(plans->pl_inv2);
     fftw_destroy_plan(plans->plan2);
     */
     fftw_cleanup();

} // end of cleanup & memory free
