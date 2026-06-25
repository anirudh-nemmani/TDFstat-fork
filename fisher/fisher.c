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
#include <time.h>
#include <dirent.h>

#include "auxi.h"
#include "settings.h"
#include "fisher.h"
#include "struct.h"
#include "init.h"

/* Fisher matrix calculation */

void fisher_matrix( Fisher_settings *sett,
                    Command_line_opts *opts,
                    Aux_arrays *aux_arr,
                    Signal_params *sig_params)
{
     int i, j, k, l;
     
     double delta_param = 1e-6; // small change in parameters for numerical derivative

     double dh_signal[7][sett->nifo][sett->N]; // derivative of the signal with respect to parameters
     double fisher_matrix[7][7]; // Fisher matrix
     double inv_fisher_matrix[7][7]; // Inverse fisher matrix

     FILE *out;

     
     // Loop over parameters
     for (i = 0; i < 7; i++) {

          Detector_settings det_pdp[MAX_DETECTORS], det_ndp[MAX_DETECTORS];
          Signal_params sig_pdp, sig_ndp;
          double dp;

          ParamIndex param = (ParamIndex)i;

          // Reset the detector settings
          if (detector_array_deep_copy(det_pdp, ifo, sett) != 0) exit(EXIT_FAILURE);
          if (detector_array_deep_copy(det_ndp, ifo, sett) != 0) exit(EXIT_FAILURE);

          
          sig_pdp = *sig_params;
          sig_ndp = *sig_params;

          switch (i) {

               case FREQ:
                    dp = sig_params->freq*delta_param;
                    sig_pdp.freq = sig_params->freq + dp;
                    sig_ndp.freq = sig_params->freq - dp;
                    break;

               case FDOT:
                    dp = sig_params->fdot*delta_param;
                    sig_pdp.fdot = sig_params->fdot + dp;
                    sig_ndp.fdot = sig_params->fdot - dp;
                    break;

               case RA:
                    dp = sig_params->ra*delta_param;
                    sig_pdp.ra = sig_params->ra + dp;
                    sig_ndp.ra = sig_params->ra - dp;
                    break;

               case DEC:
                    dp = sig_params->dec*delta_param;
                    sig_pdp.dec = sig_params->dec + dp;
                    sig_ndp.dec = sig_params->dec - dp;
                    break;

               case IOTA:
                    dp = sig_params->iota*delta_param;
                    sig_pdp.iota = sig_params->iota + dp;
                    sig_ndp.iota = sig_params->iota - dp;
                    break;

               case PSI:
                    dp = sig_params->psi*delta_param;
                    sig_pdp.psi = sig_params->psi + dp;
                    sig_ndp.psi = sig_params->psi - dp;
                    break;

               case PHASE:
                    dp = sig_params->phase*delta_param;
                    sig_pdp.phase = sig_params->phase + dp;
                    sig_ndp.phase = sig_params->phase - dp;
                    break;

               default :
                    fprintf(stderr, "Error: Invalid parameter index %d\n", i);
                    exit(EXIT_FAILURE);
          
          } // end switch

          signal_gen(sett, opts, aux_arr, &sig_pdp, det_pdp);
          signal_gen(sett, opts, aux_arr, &sig_ndp, det_ndp);

          // Compute the numerical derivative of the signal with respect to the parameter
          // Loop over detector
          for (j = 0; j < sett->nifo; j++) {
               for (k = 0; k < sett->N; k++) {
                    dh_signal[i][j][k] = (det_pdp[j].sig.signal[k] - det_ndp[j].sig.signal[k]);
               } 
          } // end loop over detector

          // Clean up
          detector_array_free_owned(det_pdp, sett);
          detector_array_free_owned(det_ndp, sett);
               
     } // loop over parameters 

     // Compute the Fisher matrix
     for (i = 0; i < 7; i++) {
          for (j = 0; j <= i; j++) {
               fisher_matrix[i][j] = 0.0;
               for (k = 0; k < sett->nifo; k++) {

                    // Trapezoidal integral method
                    double integral = 0;
                    for (l = 0; l < sett->N; l++) {
                         double w = (l == 0 || l == sett->N - 1) ? 0.5 : 1.0;
                         integral += w * dh_signal[i][k][l] * dh_signal[j][k][l];
                    }
                    fisher_matrix[i][j] += (integral * sett->dt) / (ifo[k].sig.var);
               }
               fisher_matrix[j][i] = fisher_matrix[i][j]; // Symmetric matrix
          }
     } // end Fisher loop

     // Print the Fisher matrix and save it
     char filename[562];
     sprintf(filename, "%s/fisher_matrix_%03d_%04d%s.txt", 
          opts->outdir, opts->seg, opts->band, opts->label);
     out = fopen(filename, "w");
     if (out == NULL) {
          perror(filename);
          exit(EXIT_FAILURE);
     }
     printf("\n\n\nFisher matrix:\n");
     for (i = 0; i < 7; i++) {
          for (j = 0; j < 7; j++) {
               printf("%.16e%s", fisher_matrix[i][j], (j == 6) ? "\n" : " ");
               fprintf(out, "%.16e%s", fisher_matrix[i][j], (j == 6) ? "\n" : " ");
          }
     }

     fclose(out);

     // Inverting the Fisher matrix to compute the covariance matrix
     invm(*fisher_matrix, 7, *inv_fisher_matrix);
     sprintf(filename, "%s/inv_fisher_matrix_%03d_%04d%s.txt", 
          opts->outdir, opts->seg, opts->band, opts->label);
     out = fopen(filename, "w");
     if (out == NULL) {
         perror(filename);
         exit(EXIT_FAILURE);
     }
     printf("\n\n\nInverse Fisher matrix:\n");
     for (i = 0; i < 7; i++) {
          for (j = 0; j < 7; j++) {
               printf("%.16e%s", inv_fisher_matrix[i][j], (j == 6) ? "\n" : " ");
               fprintf(out, "%.16e%s", inv_fisher_matrix[i][j], (j == 6) ? "\n" : " ");
          }
     }
     fclose(out);

} // end Fisher code


// Copy buffer function
static int copy_buf( void **dst,
                     const void *src,
                     size_t elem_size,
                     size_t count)
{
    if (!src || count == 0) {
        *dst = NULL;
        return 0;
    }

    void *p = malloc(elem_size * count);
    if (!p) return -1;

    memcpy(p, src, elem_size * count);
    *dst = p;
    return 0;
}

// Free the signal struct
static void free_signals_owned( Signals *s )
{
    free(s->noise);
    free(s->signal);
    free(s->DetSSB);
    free(s->aa);
    free(s->bb);
    free(s->shftf);
    free(s->shft);
    free(s->xDatma);
    free(s->xDatmb);
}

int detector_settings_deep_copy( Detector_settings *dst,
                                 const Detector_settings *src,
                                 Fisher_settings *sett)
{
    if (!dst || !src) return -1;

    // Start with shallow copy
    memset(dst, 0, sizeof(*dst));
    *dst = *src;

    // Copy values of 
    if (copy_buf((void**)&dst->sig.noise,  src->sig.noise,  sizeof(float),          sett->N)) goto fail;
    if (copy_buf((void**)&dst->sig.signal, src->sig.signal, sizeof(float),          sett->N)) goto fail;
    if (copy_buf((void**)&dst->sig.DetSSB, src->sig.DetSSB, sizeof(double),         sett->N)) goto fail;
    if (copy_buf((void**)&dst->sig.aa,     src->sig.aa,     sizeof(double),         sett->N)) goto fail;
    if (copy_buf((void**)&dst->sig.bb,     src->sig.bb,     sizeof(double),         sett->N)) goto fail;
    if (copy_buf((void**)&dst->sig.shftf,  src->sig.shftf,  sizeof(double),         sett->N)) goto fail;
    if (copy_buf((void**)&dst->sig.shft,   src->sig.shft,   sizeof(double),         sett->N)) goto fail;
    if (copy_buf((void**)&dst->sig.xDatma, src->sig.xDatma, sizeof(double complex), sett->N)) goto fail;
    if (copy_buf((void**)&dst->sig.xDatmb, src->sig.xDatmb, sizeof(double complex), sett->N)) goto fail;

    return 0;

fail:
    free_signals_owned(&dst->sig);
    return -1;
}

int detector_array_deep_copy( Detector_settings *dst_arr,
                              const Detector_settings *src_arr,
                              Fisher_settings *sett)
{
    if (!dst_arr || !src_arr) return -1;

    for (size_t j = 0; j < sett->nifo; ++j) {
        if (detector_settings_deep_copy(&dst_arr[j], &src_arr[j], sett) != 0) {
            // cleanup previous copied entries
            for (size_t k = 0; k < j; ++k) {
                free_signals_owned(&dst_arr[k].sig);
            }
            return -1;
        }
    }
    return 0;
}

void detector_array_free_owned( Detector_settings *arr,
                                Fisher_settings *sett)
{
    if (!arr) return;
    for (size_t j = 0; j < sett->nifo; ++j) {
        free_signals_owned(&arr[j].sig);
    }
}

void perturb_params( Signal_params *out,
                     Signal_params *base,
                     ParamIndex p,
                     double dp)
{
    *out = *base;

    switch (p) {
        case FREQ:  out->freq  += dp; break;
        case FDOT:  out->fdot  += dp; break;
        case RA:    out->ra    += dp; break;
        case DEC:   out->dec   += dp; break;
        case IOTA:  out->iota  += dp; break;
        case PSI:   out->psi   += dp; break;
        case PHASE: out->phase += dp; break;
    }
}
