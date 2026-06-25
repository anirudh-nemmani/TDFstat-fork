#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <complex.h>
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

#include "auxi.h"
#include "struct.h"
#include "init.h"
#include "settings.h"
#include "../utils/iniparser/src/iniparser.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

void read_ini_file( Fisher_settings *sett,
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

     // dictionary containing input data
     opts->indir = iniparser_getstring(ini, "fisher:indir", NULL);
     // output directory
     opts->outdir = iniparser_getstring(ini, "fisher:outdir", ".");
     // band number
     opts->band = iniparser_getint(ini, "fisher:band", 0);
     // time segment number
     opts->seg = iniparser_getint(ini, "fisher:seg", 0);
     // number of days (integer)
     sett->nod = iniparser_getint(ini, "fisher:nod", 0);
     // sampling interval of the input time series [seconds] (double)
     sett->dt = iniparser_getdouble(ini, "fisher:dt", -1.);
     // bands overlap [0-1.] defines band base frequency, fpo
     opts->overlap = iniparser_getdouble(ini, "fisher:overlap", 0.);
     // name of the file with signall to add
     opts->addsig = iniparser_getstring(ini, "fisher:addsig", "");
     // use data from subset of detectors only (default is to use all available)
     opts->usedet = iniparser_getstring(ini, "fisher:usedet", "");
     // optional label of output files
     opts->label = iniparser_getstring(ini, "fisher:label", "");
     if (strlen(opts->label)) {
         const char *tmp = opts->label;
         opts->label = malloc(strlen(tmp) + 2);
         sprintf((char*)opts->label, "_%s", tmp);
     }
     // runtime modifiers, supported values are: {read_O3}
     opts->mods = iniparser_getstring(ini, "fisher:mods", "");

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
     if (error) {
          exit(EXIT_FAILURE);
     }

     sett->fpo = 10. + (1. - opts->overlap)*opts->band*(0.5/sett->dt);

} // read_ini_file

void init_arrays( Fisher_settings *sett,
                  Command_line_opts *opts,
                  Aux_arrays *aux_arr )
{

     int i;
     size_t status;

     // Allocates and initializes to zero the data, detector ephemeris
     // and the F-statistic arrays

     FILE *data;

     for (i=0; i<sett->nifo; i++) {

          ifo[i].sig.noise = (float *) calloc(sett->N, sizeof(float));
          ifo[i].sig.signal = (float *) calloc(sett->N, sizeof(float));

          // Input time-domain data handling
          //
          // The file name ifo[i].xdatname is constructed
          // in settings.c, while looking for the detector
          // subdirectories

          if((data = fopen(ifo[i].xdatname, "r")) != NULL) {
               status = fread((void *)(ifo[i].sig.noise), sizeof(float), sett->N, data);
               fclose (data);
          } else {
               perror (ifo[i].xdatname);
               exit(EXIT_FAILURE);
          }

          int j, Nzeros=0;
          // Checking for null values in the data
          for (j=0; j < sett->N; j++)
               if(!ifo[i].sig.noise[j]) Nzeros++;

          ifo[i].sig.Nzeros = Nzeros;

          // factor N/(N - Nzeros) to account for null values in the data
          ifo[i].sig.crf0 = (double)sett->N/(sett->N - ifo[i].sig.Nzeros);

          // Estimation of the variance for each detector
          ifo[i].sig.var = (ifo[i].sig.crf0)*var(ifo[i].sig.noise, sett->N);

          ifo[i].sig.DetSSB = (double *) calloc(3*sett->N, sizeof(double));

          // Ephemeris file handling
          char filename[562];
          sprintf (filename, "%s/%03d/%s/DetSSB.bin", opts->indir, opts->seg, ifo[i].name);

          if((data = fopen(filename, "r")) != NULL) {
               // Detector position w.r.t Solar System Baricenter
               // for every datapoint
               status = fread((void *)(ifo[i].sig.DetSSB), sizeof(double), 3*sett->N, data);

               // Deterministic phase defining the position of the Earth
               // in its diurnal motion at t=0
               status = fread((void *)(&ifo[i].sig.phir), sizeof(double), 1, data);

               // Earth's axis inclination to the ecliptic at t=0
               status = fread((void *)(&ifo[i].sig.epsm), sizeof(double), 1, data);
               fclose (data);

          } else {
               perror (filename);
               return ;
          }

          // Start time reading
          sprintf (filename, "%s/%03d/%s/starting_date", opts->indir, opts->seg, ifo[i].name);

          if ((data = fopen(filename, "r")) != NULL) {
               // Start time of the data segment in GPS seconds
               status = fscanf(data, "%lf", &ifo[i].start_time);
               fclose (data);
               printf("[%s] Starting time = %.3f\n", ifo[i].name, ifo[i].start_time);
          } else {
               perror (filename);
               return ;
          }

          // sincos
          ifo[i].sig.sphir = sin(ifo[i].sig.phir);
          ifo[i].sig.cphir = cos(ifo[i].sig.phir);
          ifo[i].sig.sepsm = sin(ifo[i].sig.epsm);
          ifo[i].sig.cepsm = cos(ifo[i].sig.epsm);

          sett->sepsm = ifo[i].sig.sepsm;
          sett->cepsm = ifo[i].sig.cepsm;

          ifo[i].sig.aa = (double *) calloc(sett->N, sizeof(double));
          ifo[i].sig.bb = (double *) calloc(sett->N, sizeof(double));

          ifo[i].sig.shft = (double *) calloc(sett->N, sizeof(double));
          ifo[i].sig.shftf = (double *) calloc(sett->N, sizeof(double));

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

     // Auxiliary arrays, Earth's rotation
     aux_arr->t2 = (double *) calloc(sett->N, sizeof (double));
     aux_arr->cosmodf = (double *) calloc(sett->N, sizeof (double));
     aux_arr->sinmodf = (double *) calloc(sett->N, sizeof (double));
     double omrt;

     for (i=0; i<sett->N; i++) {
          omrt = (sett->omr)*i;     // Earth angular velocity * dt * i
          aux_arr->t2[i] = sqr((double)i);
          aux_arr->cosmodf[i] = cos(omrt);
          aux_arr->sinmodf[i] = sin(omrt);
     }

} // end of init arrays

/* Read signal parameters from a file */

void read_signal_file( Signal_params *sgnl_params,
                       Command_line_opts *opts)
{
     FILE *data;
     char amporsnr[4];
     
     if ((data=fopen (opts->addsig, "r")) != NULL) {
         // Fscanning for the GW amplitude h0 or signal-to-noise,
         // the grid size and the reference frame
         // (for which the signal freq. is not spun-down/up)
         
         do {
             fscanf (data, "%s", amporsnr);
         } while ( strcmp(amporsnr, "amp")!=0 && strcmp(amporsnr, "snr")!=0 );
         
         strcpy(sgnl_params->amporsnr, amporsnr);
         if(!strcmp(amporsnr, "amp")) {
             fscanf (data, "%le %d %le %le %le %le %le %le %le",
                                  &sgnl_params->h0, &sgnl_params->reffr,
                                  &sgnl_params->freq, &sgnl_params->fdot, &sgnl_params->ra, &sgnl_params->dec,
                                  &sgnl_params->iota, &sgnl_params->psi, &sgnl_params->phase);
             
             printf("add_signal(): GW amplitude h0 is %le\n   The reference band of the signal is %d\n"
                    "   The signal is injected at the following parameters:\n"
                    "   Frequency [Hz]         : %le\n"
                    "   Spin-down [Hz/s]       : %le\n"
                    "   Right ascension [rad]  : %le\n"
                    "   Declination [rad]      : %le\n"
                    "   Inclination [rad]      : %le\n"
                    "   Polarization [rad]     : %le\n"
                    "   Phase [rad]            : %le\n",
                    sgnl_params->h0, sgnl_params->reffr, sgnl_params->freq, sgnl_params->fdot, sgnl_params->ra, sgnl_params->dec,
                    sgnl_params->iota, sgnl_params->psi, sgnl_params->phase);
         } else if (!strcmp(amporsnr, "snr")) {
             fscanf (data, "%le %d %le %le %le %le %le %le %le",
                                  &sgnl_params->snr, &sgnl_params->reffr,
                                  &sgnl_params->freq, &sgnl_params->fdot, &sgnl_params->ra, &sgnl_params->dec,
                                  &sgnl_params->iota, &sgnl_params->psi, &sgnl_params->phase);
             
             printf("add_signal(): GW (network) signal-to-noise ratio is %le\n   The reference band of the signal is %d\n"
                    "   The signal is injected at the following parameters:\n"
                    "   Frequency [Hz]         : %le\n"
                    "   Spin-down [Hz/s]       : %le\n"
                    "   Right ascension [rad]  : %le\n"
                    "   Declination [rad]      : %le\n"
                    "   Inclination [rad]      : %le\n"
                    "   Polarization [rad]     : %le\n"
                    "   Phase [rad]            : %le\n",
                    sgnl_params->snr, sgnl_params->reffr, sgnl_params->freq, sgnl_params->fdot, sgnl_params->ra, sgnl_params->dec,
                    sgnl_params->iota, sgnl_params->psi, sgnl_params->phase);
         } else {
             printf("Invalid format in signal file. First column of signals should start with 'amp' or 'snr'.\n");
             exit(0);
         }
         fclose (data);
         
     } else {
         perror (opts->addsig);
     }
} // end of read_signal_file



/* Add signal to data */

void signal_gen( Fisher_settings *sett,
                 Command_line_opts *opts,
                 Aux_arrays *aux_arr,
                 Signal_params *sgnl_params,
                 Detector_settings *ifo)
{

     int i, j, n, gsize, reffr;
     double sum = 0., cof, d1;
     double sigma_noise = 1.0;
     double be[2];
     double sinalpha, cosalpha, sindelta, cosdelta, phaseadd, shiftadd;
     double phi, psi, cosi, cosip, iota, amplit[4];
     double nSource[3], freqo[2];

     char amporsnr[4];
          
     // Setting the reference frame
     reffr = sgnl_params->reffr;
     
     // Converting the frequency into dimensionless units
     freqo[0] = 2 * M_PI * sgnl_params->freq * sett->dt;
     freqo[1] = M_PI * sgnl_params->fdot * sett->dt * sett->dt;

     // Shift the frequency to the current segment based on spindown from the reference segment
     freqo[0] += 2.*freqo[1]*(sett->N)*(opts->seg - reffr);

     // Check if the signal is in band
     if( freqo[0]/(2*M_PI*sett->dt) < sett->fpo || freqo[0]/(2*M_PI*sett->dt) > sett->fpo + sett->B ) {
          printf("add_signal(): signal out of band f=%le s=%le\n", freqo[0]/(2*M_PI*sett->dt), sgnl_params->fdot);
          return;
     }

     // Calculation of sin alpha, cos alpha, sin delta, cos delta of the signal.
     // Check Eq. 18 of Phys. Rev. D 58, 063001 1998
     sinalpha = sin(sgnl_params->ra);
     cosalpha = cos(sgnl_params->ra);
     sindelta = sin(sgnl_params->dec);
     cosdelta = cos(sgnl_params->dec);

     // Calculation of four amplitudes from polarization, phase and inclination
     // Check Eq. 32 - 35 of Phys. Rev. D 58, 063001 1998
     cosi = cos(sgnl_params->iota);
     psi = sgnl_params->psi;
     phi = sgnl_params->phase;
     cosip = (1. + cosi*cosi)/2.;

     amplit[0] = cos(2.*psi)*cosip*cos(phi) - sin(2.*psi)*cosi*sin(phi);
     amplit[1] = sin(2.*psi)*cosip*cos(phi) + cos(2.*psi)*cosi*sin(phi);
     amplit[2] = -cos(2.*psi)*cosip*sin(phi) - sin(2.*psi)*cosi*cos(phi);
     amplit[3] = -sin(2.*psi)*cosip*sin(phi) + cos(2.*psi)*cosi*cos(phi);

     // To keep coherent phase between time segments
     double phaseshift = freqo[0]*sett->N*(opts->seg - reffr)
                       - freqo[1]*pow(sett->N*(opts->seg - reffr), 2);

     // Allocate arrays for added signal, for each detector
     double **signadd = malloc((sett->nifo)*sizeof(double *));
     for (n=0; n<sett->nifo; n++)
          signadd[n] = malloc((sett->N)*sizeof(double));

     // Loop for each detector - sum calculations
     for (n=0; n<sett->nifo; n++) {

          modvir(sinalpha, cosalpha, sindelta, cosdelta, sett->N, &ifo[n], aux_arr);

          nSource[0] = cosalpha*cosdelta;
          nSource[1] = sinalpha*cosdelta;
          nSource[2] = sindelta;

          for (i=0; i<sett->N; i++) {

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
               signadd[n][i] = amplit[0]*(ifo[n].sig.aa[i])*cos(phaseadd)
                             + amplit[1]*(ifo[n].sig.bb[i])*cos(phaseadd)
                             + amplit[2]*(ifo[n].sig.aa[i])*sin(phaseadd)
                             + amplit[3]*(ifo[n].sig.bb[i])*sin(phaseadd);

               // Sum over signals
               sum += pow(signadd[n][i], 2.);

          } // data loop
     } // detector loop

     // Signal amplitude h0 from the snr
     // (currently only makes sense for Gaussian noise with fixed sigma and one detector, but can be generalized)
     if (!strcmp(sgnl_params->amporsnr, "snr")) {
         sgnl_params->h0 = (sgnl_params->snr * sqrt(ifo[0].sig.var)) / (sqrt(sum));
     }

     // Loop for each detector - adding signal to data (point by point)
     for (n=0; n<sett->nifo; n++) {
          for (i=0; i<sett->N; i++) {
               // Adding the signal to the data vector
               if (ifo[n].sig.noise[i])
                    ifo[n].sig.signal[i] = sgnl_params->h0*signadd[n][i];
          } // data loop
     } // detector loop

     // Free auxiliary 2d array
     for (n=0; n<sett->nifo; n++)
          free(signadd[n]);
     free(signadd);

} // add_signal()

/* Cleanup & memory free  */

void cleanup( Fisher_settings *sett,
            Command_line_opts *opts,
            Aux_arrays *aux)
{

   int i;

   for(i=0; i<sett->nifo; i++) {
        free(ifo[i].sig.noise);
        free(ifo[i].sig.signal);
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

} // end of cleanup & memory free
