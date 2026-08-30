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

#define NPARAMS 7

static double get_param(const Signal_params *p, ParamIndex idx)
{
     switch (idx) {
          case FREQ:  return p->freq;
          case FDOT:  return p->fdot;
          case RA:    return p->ra;
          case DEC:   return p->dec;
          case IOTA:  return p->iota;
          case PSI:   return p->psi;
          case PHASE: return p->phase;
     }
     return 0.0;
}

void perturb_params( Signal_params *out,
                     const Signal_params *base,
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

// Finite-difference step sizes.

// The step for each parameter must be small compared to that parameter's
// NATURAL scale (roughly the width of the likelihood), not a fixed relative
// fraction of its value:

//   - frequency resolution over an observation of length T = N*dt is ~1/T,
//     so dp_freq << 1/T (a relative step f*1e-6 changes the accumulated
//     phase by tens of radians over 6 days and destroys the derivative);
//   - spindown resolution is ~1/T^2;
//   - the angles enter through O(1) trigonometric factors (iota, psi,
//     phase) or through the Doppler phase ~2*pi*f*(AU/c) (ra, dec), so a
//     small absolute step in radians is appropriate.

static void choose_steps(const Fisher_settings *sett, double dp[NPARAMS])
{
     const double T = sett->N * sett->dt;   /* observation time [s] */

     dp[FREQ]  = 1.0e-3 / T;                /* [Hz]   */
     dp[FDOT]  = 1.0e-3 / (T*T);            /* [Hz/s] */
     dp[RA]    = 1.0e-6;                    /* [rad]  */
     dp[DEC]   = 1.0e-6;                    /* [rad]  */
     dp[IOTA]  = 1.0e-6;                    /* [rad]  */
     dp[PSI]   = 1.0e-6;                    /* [rad]  */
     dp[PHASE] = 1.0e-6;                    /* [rad]  */
}

static void save_matrix( const double *m, int n, const char *heading,
                         const char *fnprefix, const Command_line_opts *opts)
{
     char filename[FNAME_LENGTH + 64];
     snprintf(filename, sizeof(filename), "%s/%s_%03d_%04d%s.txt",
              opts->outdir, fnprefix, opts->seg, opts->band, opts->label);

     FILE *out = fopen(filename, "w");
     if (!out) {
          perror(filename);
          exit(EXIT_FAILURE);
     }

     printf("\n%s:\n", heading);
     for (int i = 0; i < n; i++) {
          for (int j = 0; j < n; j++) {
               const char *sep = (j == n-1) ? "\n" : " ";
               printf("%.16e%s", m[n*i + j], sep);
               fprintf(out, "%.16e%s", m[n*i + j], sep);
          }
     }
     fclose(out);
}

// Fisher matrix
// F_ij = sum_det (1/sigma_det^2) * Integral dt  (dh/dp_i)(dh/dp_j)
// Derivatives are central finite differences:
//      dh/dp_i = ( h(p_i + dp_i) - h(p_i - dp_i) ) / (2 dp_i)

// signal_gen() fully overwrites ifo[n].sig.aa, .bb and .signal on each
// call and reads nothing that a previous call modified, so it is safe to
// evaluate all perturbed signals in place on the global ifo[] array and
// snapshot the result -- no deep copies of the detector structs needed.

void fisher_matrix( Fisher_settings *sett,
                    Command_line_opts *opts,
                    Aux_arrays *aux_arr,
                    Signal_params *sig_params)
{
     const int nifo = sett->nifo;
     const int N    = sett->N;

     if (!strcmp(sig_params->amporsnr, "snr")) {
          signal_gen(sett, opts, aux_arr, sig_params, ifo); /* sets h0 */
          strcpy(sig_params->amporsnr, "amp");
     }

     double *dh = malloc((size_t)NPARAMS * nifo * N * sizeof *dh);
     float  *hp = malloc((size_t)nifo * N * sizeof *hp);
     if (!dh || !hp) {
          fprintf(stderr, "fisher_matrix: out of memory\n");
          exit(EXIT_FAILURE);
     }
#define DH(i,d,t) dh[((size_t)(i)*nifo + (d))*N + (t)]
#define HP(d,t)   hp[(size_t)(d)*N + (t)]

     double dp[NPARAMS];
     choose_steps(sett, dp);

     for (int i = 0; i < NPARAMS; i++) {

          Signal_params sig_p, sig_n;
          perturb_params(&sig_p, sig_params, (ParamIndex)i, +dp[i]);
          perturb_params(&sig_n, sig_params, (ParamIndex)i, -dp[i]);

          // h(p + dp), snapshot it
          signal_gen(sett, opts, aux_arr, &sig_p, ifo);
          for (int d = 0; d < nifo; d++)
               for (int t = 0; t < N; t++)
                    HP(d,t) = ifo[d].sig.signal[t];

          // h(p - dp), form the central difference
          signal_gen(sett, opts, aux_arr, &sig_n, ifo);
          const double inv2dp = 1.0 / (2.0 * dp[i]);
          for (int d = 0; d < nifo; d++)
               for (int t = 0; t < N; t++)
                    DH(i,d,t) = (HP(d,t) - ifo[d].sig.signal[t]) * inv2dp;

          printf("Parameter %d: step %.3e, derivative computed\n", i, dp[i]);
     }

     free(hp);

     // Leave ifo[] holding the signal at the (unperturbed) base parameters 
     signal_gen(sett, opts, aux_arr, sig_params, ifo);

     // F_ij: trapezoidal time integral, whitened per detector, summed
     // over the network.  Symmetric, so compute the lower triangle.
     double F[NPARAMS][NPARAMS], C[NPARAMS][NPARAMS];

     for (int i = 0; i < NPARAMS; i++) {
          for (int j = 0; j <= i; j++) {
               double fij = 0.0;
               for (int d = 0; d < nifo; d++) {
                    double integral = 0.0;
                    for (int t = 0; t < N; t++) {
                         double w = (t == 0 || t == N-1) ? 0.5 : 1.0;
                         integral += w * DH(i,d,t) * DH(j,d,t);
                    }
                    fij += integral * sett->dt / ifo[d].sig.var;
               }
               F[i][j] = F[j][i] = fij;
          }
     }

     free(dh);
#undef DH
#undef HP

     save_matrix(&F[0][0], NPARAMS, "Fisher matrix", "fisher_matrix", opts);

     // Covariance matrix = inverse of the Fisher matrix
     if (invm(&F[0][0], NPARAMS, &C[0][0])) {
          fprintf(stderr, "fisher_matrix: Fisher matrix is singular, "
                          "covariance matrix not computed\n");
          return;
     }
     save_matrix(&C[0][0], NPARAMS, "Inverse Fisher matrix (covariance)",
                 "inv_fisher_matrix", opts);

} // fisher_matrix
