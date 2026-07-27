// general purpose functions for the project

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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
