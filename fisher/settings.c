#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <float.h>

#include "auxi.h"
#include "struct.h"
#include "settings.h"
#include <glob.h>

/* Fisher settings:
 * bandwidth and Earth parameters and other details
 */

void fisher_settings(Fisher_settings* sett)
{

     double dt, B, oms, omr;
     int N, s, nd;

     dt = sett->dt;                    // data sampling time:
     B = 0.5/dt;                       // Bandwidth
     oms = 2.*M_PI*(sett->fpo)*dt;     // Dimensionless angular frequency
     omr = C_OMEGA_R*dt;

     N = round (sett->nod*C_SIDDAY/dt);      // No. of data points

     s = 1;                                    // No. of spindowns

     nd = 2;     // Degree of freedom, (2*nd = deg. no ofrees of freedom for chi^2)

     sett->B=B;          	// bandwidth
     sett->oms=oms;      	// dimensionless angular frequency
     sett->omr=omr;      	// C_OMEGA_R * dt
     sett->N=N;          	// number of data points
     sett->s=s;          	// number of spindowns
     sett->nd=nd;        	// degrees of freedom

     // calculate 1/day frequency in units of Fstat bins
     // signal width is ~ 5*dd
     double dayf = 1./C_SIDDAY;           // 1/day frequency
     
     printf("------------------------ Settings --------------------------\n");
     printf(" B         N         fpo         NAV\n");
     printf(" %-9.3f %-12d %-9.3f %-11d\n", sett->B, sett->N, sett->fpo, NAV);
     printf("------------------------------------------------------------\n");


} // fisher settings


/* Network of detectors' discovery:
* finds subdirectories in the main input directory,
* which by convention should be named like V1, L1, H1
* and which contain input data and ephemerids;
* writes appropriate detector-related data into structs.
*/

void detectors_settings( Fisher_settings* sett, Command_line_opts *opts)
{

     int i=0, j=0;
     char dirname[1024], x[1332];
     DIR *dp;
     FILE *data;
     const char *dets[] = {"H1", "L1", "V1"};
     char det[DETNAME_LENGTH];

     // Test frame input directory
     sprintf (dirname, "%s/%03d", opts->indir, opts->seg);
     dp = opendir(dirname);
     if (dp) {
          closedir(dp);
     } else {
          printf("Can't open the input directory: %s", dirname);
          exit(EXIT_FAILURE);
     }

     // test availability of data for detectors
     for (i=0; i<3; i++) {
          if ( !strlen(opts->usedet) || (strlen(opts->usedet) && (strstr(opts->usedet, dets[i]))) ) {
               // detector directory
               memset(dirname, 0, sizeof(dirname));
               sprintf (dirname, "%s/%03d/%s", opts->indir, opts->seg, dets[i]);
               dp = opendir(dirname);
               if (dp) {
                    closedir(dp);
                    sprintf(x, "%s/xdat_%03d_%04d.bin", dirname, opts->seg, opts->band);
                    data = fopen(x, "r");
                    if (data) {
                         fclose(data);
                         //strncpy(ifo[j].xdatname, x, strlen(x));
                         strcpy(ifo[j].xdatname, x);
                         strncpy(ifo[j].name, dets[i], DETNAME_LENGTH);
                         memset(x, 0, sizeof(x));
                         j++;
                    } else {
                         printf("Directory %s exists, but no input file found:\n%s missing...\n", dirname, x);
                         exit(EXIT_FAILURE);
                    }
               } else {
                    if ( strlen(opts->usedet) && (strstr(opts->usedet, dets[i])) ) {
                         printf("Can't open the input directory requied by -usedet: %s", dirname);
                         exit(EXIT_FAILURE);
                    }
               } // if dp
          } // if
     } // for i

     sett->nifo = j;

     for(i=0; i<sett->nifo; i++) {

          //    printf("Using %s IFO as detector #%d... %s as input time series data\n",
          printf("IFO[%d] = %s , data = %s\n", i, ifo[i].name, ifo[i].xdatname);

          if(!strcmp("V1", ifo[i].name)) {
               // Virgo detector

               // Geographical latitude phi in radians
               ifo[i].ephi = (43.+37./60.+53.0880/3600.)/RAD_TO_DEG;
               // Geographical longitude in radians
               ifo[i].elam = (10.+30./60.+16.1885/3600.)/RAD_TO_DEG;
               // Height h above the Earth ellipsoid in meters
               ifo[i].eheight = 51.884;
               // Orientation of the detector gamma
               ifo[i].egam = (135. - (19.0+25./60.0+57.96/3600.))/RAD_TO_DEG;

          } else if(!strcmp("H1", ifo[i].name )) {
               // Hanford H1 detector

               // Geographical latitude phi in radians
               ifo[i].ephi = (46+(27+18.528/60.)/60.)/RAD_TO_DEG;
               // Geographical longitude in radians
               ifo[i].elam = -(119+(24+27.5657/60.)/60.)/RAD_TO_DEG;
               // Height h above the Earth ellipsoid in meters
               ifo[i].eheight = 142.554;
               // Orientation of the detector gamma
               ifo[i].egam = 170.9994/RAD_TO_DEG;

          } else if(!strcmp("L1", ifo[i].name )) {
               // Livingston L1 detector

               // Geographical latitude phi in radians
               ifo[i].ephi = (30+(33+46.4196/60.)/60.)/RAD_TO_DEG;
               // Geographical longitude in radians
               ifo[i].elam = -(90+(46+27.2654/60.)/60.)/RAD_TO_DEG;
               // Height h above the Earth ellipsoid in meters
               ifo[i].eheight = -6.574;
               // Orientation of the detector gamma
               ifo[i].egam = 242.7165/RAD_TO_DEG;
          }
     }  // for i

     // todo: check if there are -usedet detectors without directory match

} // detectors settings

/* Coefficients of the amplitude modulation functions
* of the Virgo detector
*/

void rogcvir(Detector_settings *ifo)
{

     /* In the notation of Phys. Rev. D 58, 063001 (1998):
     * ephi = lambda (geographical latitude phi in radians)
     * egam = gamma (orientation of the detector)
     *
     * (see modvir function in jobcore.c for Eqs. 12 and 13)
     */

     //printf("Calculating the amplitude modulation functions for %s...\n", ifo->name);

     ifo->amod.c1 = .25*sin(2.*ifo->egam)*(1+sqr(sin(ifo->ephi)));
     ifo->amod.c2 = -.5*cos(2.*ifo->egam)*sin(ifo->ephi);
     ifo->amod.c3 = .5*sin(2.*ifo->egam)*sin(2.*ifo->ephi);
     ifo->amod.c4 = -cos(2.*ifo->egam)*cos(ifo->ephi);
     ifo->amod.c5 = .75*sin(2.*ifo->egam)*sqr(cos(ifo->ephi));
     ifo->amod.c6 = cos(2.*ifo->egam)*sin(ifo->ephi);
     ifo->amod.c7 = .5*sin(2.*ifo->egam)*(1.+sqr(sin(ifo->ephi)));
     ifo->amod.c8 = cos(2.*ifo->egam)*cos(ifo->ephi);
     ifo->amod.c9 = .5*sin(2.*ifo->egam)*sin(2.*ifo->ephi);

} // rogcvir


/* Amplitude modulation of the signal
*/

void modvir( double sinal, double cosal, double sindel, double cosdel,
             int Np, Detector_settings *ifo, Aux_arrays *aux)
{

     int t;
     double cosalfr, sinalfr, c2d, c2sd, c, s, c2s, cs;

     double c1 = ifo->amod.c1,
            c2 = ifo->amod.c2,
            c3 = ifo->amod.c3,
            c4 = ifo->amod.c4,
            c5 = ifo->amod.c5,
            c6 = ifo->amod.c6,
            c7 = ifo->amod.c7,
            c8 = ifo->amod.c8,
            c9 = ifo->amod.c9;

     cosalfr = cosal*(ifo->sig.cphir) + sinal*(ifo->sig.sphir);
     sinalfr = sinal*(ifo->sig.cphir) - cosal*(ifo->sig.sphir);
     c2d = sqr(cosdel);
     c2sd = sindel*cosdel;

     // Modulation factors aa, bb for every NON-ZERO data point
     for (t=0; t<Np; t++) {
          if ( fabs(ifo->sig.noise[t]) > DBL_MIN ) {
               c = cosalfr*aux->cosmodf[t] + sinalfr*aux->sinmodf[t];
               s = sinalfr*aux->cosmodf[t] - cosalfr*aux->sinmodf[t];
               c2s = 2.*sqr(c);
               cs = c*s;

               ifo->sig.aa[t] = c1*(2.-c2d)*c2s + c2*(2.-c2d)*2.*cs +
                    c3*c2sd*c + c4*c2sd*s - c1*(2.-c2d) + c5*c2d;

               ifo->sig.bb[t] = c6*sindel*c2s + c7*sindel*2.*cs +
                    c8*cosdel*c + c9*cosdel*s - c6*sindel;
          } else {
               ifo->sig.aa[t] = 0.;
               ifo->sig.bb[t] = 0.;
          }
     }
} // modvir
