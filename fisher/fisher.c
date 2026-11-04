/*
    This code contains the fisher matrix calculations for the 2f component model in TDFstat signal generation.

    The derivatives are done analytically in mathematica and then later converted into C.

    To compute the inverse of the Fisher matrix, we invert the normalised form of the Fisher matrix, R
    R = D^-1 F D^-1, D = diag(sqrt(F_jj))

    The covariance matrix, C = F^-1 = D^-1 R^-1 D^-1

    Parameters (ParamIndex order): freq, fdot, ra, dec, iota, psi, phase, h0.

    This script also contains the systematics from Cutler and Vallisneri 2018 (https://arxiv.org/pdf/0707.2982)
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
 
#include "auxi.h"
#include "settings.h"
#include "fisher.h"
#include "struct.h"
#include "init.h"
 
typedef struct {
     // signal parameters
     double alpha, delta, f, fdot;
     // amplitudes and their derivatives
     double A21, A22, A23, A24;
     double dA21_dh0,   dA22_dh0,   dA23_dh0,   dA24_dh0;    // amplit[0..3]
     double dA21_diota, dA22_diota, dA23_diota, dA24_diota;
     double dA21_dpsi,  dA22_dpsi,  dA23_dpsi,  dA24_dpsi;
     double dA21_dphi0, dA22_dphi0, dA23_dphi0, dA24_dphi0;
     // phase constants
     double dt, K, freqo0, freqo1, phaseshift, wh, omr;
     // detector: rogcvir() coefficients and phase phi_r
     double c1, c2, c3, c4, c5, c6, c7, c8, c9, phir;
} Deriv_ctx;
 
 
// Amplitudes A2i (JKS Eqs. 32-35, theta = zeta = pi/2) and their derivatives.
// They depend only on h0, iota, psi, phi0, so they are computed once.
static void amplitude_derivs(Deriv_ctx *c, double h0, double iota, double psi, double phi0)
{
     // 1. amplitudes
     c->A21 = h0*(0.5*(1 + pow(cos(iota),2))*cos(phi0)*cos(2*psi) - cos(iota)*sin(phi0)*sin(2*psi));
     c->A22 = h0*(cos(iota)*cos(2*psi)*sin(phi0) + 0.5*(1 + pow(cos(iota),2))*cos(phi0)*sin(2*psi));
     c->A23 = h0*(-0.5*(1 + pow(cos(iota),2))*cos(2*psi)*sin(phi0) - cos(iota)*cos(phi0)*sin(2*psi));
     c->A24 = h0*(cos(iota)*cos(phi0)*cos(2*psi) - 0.5*(1 + pow(cos(iota),2))*sin(phi0)*sin(2*psi));
 
     // amplitude derivatives
     c->dA21_dh0 = 0.5*(1 + pow(cos(iota),2))*cos(phi0)*cos(2*psi) - cos(iota)*sin(phi0)*sin(2*psi);
     c->dA22_dh0 = cos(iota)*cos(2*psi)*sin(phi0) + 0.5*(1 + pow(cos(iota),2))*cos(phi0)*sin(2*psi);
     c->dA23_dh0 = -0.5*(1 + pow(cos(iota),2))*cos(2*psi)*sin(phi0) - cos(iota)*cos(phi0)*sin(2*psi);
     c->dA24_dh0 = cos(iota)*cos(phi0)*cos(2*psi) - 0.5*(1 + pow(cos(iota),2))*sin(phi0)*sin(2*psi);
 
     c->dA21_diota = h0*(-(cos(iota)*cos(phi0)*cos(2*psi)*sin(iota)) + sin(iota)*sin(phi0)*sin(2*psi));
     c->dA22_diota = h0*(-(cos(2*psi)*sin(iota)*sin(phi0)) - cos(iota)*cos(phi0)*sin(iota)*sin(2*psi));
     c->dA23_diota = h0*(cos(iota)*cos(2*psi)*sin(iota)*sin(phi0) + cos(phi0)*sin(iota)*sin(2*psi));
     c->dA24_diota = h0*(-(cos(phi0)*cos(2*psi)*sin(iota)) + cos(iota)*sin(iota)*sin(phi0)*sin(2*psi));
 
     c->dA21_dpsi = h0*(-2*cos(iota)*cos(2*psi)*sin(phi0) - (1 + pow(cos(iota),2))*cos(phi0)*sin(2*psi));
     c->dA22_dpsi = h0*((1 + pow(cos(iota),2))*cos(phi0)*cos(2*psi) - 2*cos(iota)*sin(phi0)*sin(2*psi));
     c->dA23_dpsi = h0*(-2*cos(iota)*cos(phi0)*cos(2*psi) + (1 + pow(cos(iota),2))*sin(phi0)*sin(2*psi));
     c->dA24_dpsi = h0*(-((1 + pow(cos(iota),2))*cos(2*psi)*sin(phi0)) - 2*cos(iota)*cos(phi0)*sin(2*psi));
 
     c->dA21_dphi0 = h0*(-0.5*(1 + pow(cos(iota),2))*cos(2*psi)*sin(phi0) - cos(iota)*cos(phi0)*sin(2*psi));
     c->dA22_dphi0 = h0*(cos(iota)*cos(phi0)*cos(2*psi) - 0.5*(1 + pow(cos(iota),2))*sin(phi0)*sin(2*psi));
     c->dA23_dphi0 = h0*(-0.5*(1 + pow(cos(iota),2))*cos(phi0)*cos(2*psi) + cos(iota)*sin(phi0)*sin(2*psi));
     c->dA24_dphi0 = h0*(-(cos(iota)*cos(2*psi)*sin(phi0)) - 0.5*(1 + pow(cos(iota),2))*cos(phi0)*sin(2*psi));
}
 
 
// Signal value h and its 7 derivatives at sample i.
// r = &DetSSB[3*i]  (rx, ry, rz in units of c*dt).
static inline void signal_derivs(const Deriv_ctx *c, int i, const double *r,
                                 double *h, double dh[NPARAMS])
{
     const double alpha = c->alpha, delta = c->delta;
     const double rx = r[0], ry = r[1], rz = r[2];
 
     // 4. Doppler term n.r (shiftadd) and its derivatives
     double nr        = rx*cos(alpha)*cos(delta) + ry*cos(delta)*sin(alpha) + rz*sin(delta);
     double dnrDalpha = ry*cos(alpha)*cos(delta) - rx*cos(delta)*sin(alpha);
     double dnrDdelta = rz*cos(delta) - rx*cos(alpha)*sin(delta) - ry*sin(alpha)*sin(delta);
 
     // 3. modulation functions, exactly as in modvir():
     //    cc = cos x, ss = sin x, x = alpha - phi_r - omr*i
     const double x    = alpha - c->phir - c->omr*i;
     const double cc   = cos(x), ss = sin(x);
     const double sd   = sin(delta), cd = cos(delta);
     const double c2d  = cd*cd, c2sd = sd*cd, c2del = cos(2*delta);
 
     double a = c->c1*(2 - c2d)*2*cc*cc + c->c2*(2 - c2d)*2*cc*ss + c->c3*c2sd*cc + c->c4*c2sd*ss
              - c->c1*(2 - c2d) + c->c5*c2d;
     double b = c->c6*sd*2*cc*cc + c->c7*sd*2*cc*ss + c->c8*cd*cc + c->c9*cd*ss - c->c6*sd;
 
     double daDalpha = -4*c->c1*(2 - c2d)*cc*ss + 2*c->c2*(2 - c2d)*(cc*cc - ss*ss)
                     - c->c3*c2sd*ss + c->c4*c2sd*cc;
     double daDdelta = 4*c->c1*c2sd*cc*cc + 4*c->c2*c2sd*cc*ss + c->c3*c2del*cc + c->c4*c2del*ss
                     - 2*c->c1*c2sd - 2*c->c5*c2sd;
     double dbDalpha = -4*c->c6*sd*cc*ss + 2*c->c7*sd*(cc*cc - ss*ss) - c->c8*cd*ss + c->c9*cd*cc;
     double dbDdelta = 2*c->c6*cd*cc*cc + 2*c->c7*cd*cc*ss - c->c8*sd*cc - c->c9*sd*ss - c->c6*cd;
 
     // phase, exactly as phaseadd in signal_gen(), and its derivatives
     double Phi = c->freqo0*(i + nr) + c->freqo1*(2.*nr*i + (double)i*i + nr*nr)
                + c->phaseshift - c->wh*i;
     double cosPhi = cos(Phi), sinPhi = sin(Phi);
 
     const double tK = i + nr + c->K;                     // t_SSB / dt from the reference segment
     double dPhi_df    = 2*M_PI*c->dt*tK;
     double dPhi_dfdot = M_PI*c->dt*c->dt*tK*tK;
     double dPhi_dnr   = 2*M_PI*c->dt*(c->f + c->fdot*c->dt*tK);
 
     // cw modes
     double h21 = a*cosPhi;
     double h22 = b*cosPhi;
     double h23 = a*sinPhi;
     double h24 = b*sinPhi;
 
     // 2. mode derivatives
     double dh21_df = -(a*sinPhi*dPhi_df);
     double dh22_df = -(b*sinPhi*dPhi_df);
     double dh23_df =   a*cosPhi*dPhi_df;
     double dh24_df =   b*cosPhi*dPhi_df;
 
     double dh21_dfdot = -(a*sinPhi*dPhi_dfdot);
     double dh22_dfdot = -(b*sinPhi*dPhi_dfdot);
     double dh23_dfdot =   a*cosPhi*dPhi_dfdot;
     double dh24_dfdot =   b*cosPhi*dPhi_dfdot;
 
     double dh21_dalpha = cosPhi*daDalpha - a*sinPhi*dPhi_dnr*dnrDalpha;
     double dh22_dalpha = cosPhi*dbDalpha - b*sinPhi*dPhi_dnr*dnrDalpha;
     double dh23_dalpha = sinPhi*daDalpha + a*cosPhi*dPhi_dnr*dnrDalpha;
     double dh24_dalpha = sinPhi*dbDalpha + b*cosPhi*dPhi_dnr*dnrDalpha;
 
     double dh21_ddelta = cosPhi*daDdelta - a*sinPhi*dPhi_dnr*dnrDdelta;
     double dh22_ddelta = cosPhi*dbDdelta - b*sinPhi*dPhi_dnr*dnrDdelta;
     double dh23_ddelta = sinPhi*daDdelta + a*cosPhi*dPhi_dnr*dnrDdelta;
     double dh24_ddelta = sinPhi*dbDdelta + b*cosPhi*dPhi_dnr*dnrDdelta;
 
     // signal and full derivatives of h
     *h = c->A21*h21 + c->A22*h22 + c->A23*h23 + c->A24*h24;
 
     dh[FREQ]  = c->A21*dh21_df     + c->A22*dh22_df     + c->A23*dh23_df     + c->A24*dh24_df;
     dh[FDOT]  = c->A21*dh21_dfdot  + c->A22*dh22_dfdot  + c->A23*dh23_dfdot  + c->A24*dh24_dfdot;
     dh[RA]    = c->A21*dh21_dalpha + c->A22*dh22_dalpha + c->A23*dh23_dalpha + c->A24*dh24_dalpha;
     dh[DEC]   = c->A21*dh21_ddelta + c->A22*dh22_ddelta + c->A23*dh23_ddelta + c->A24*dh24_ddelta;
     dh[IOTA]  = c->dA21_diota*h21 + c->dA22_diota*h22 + c->dA23_diota*h23 + c->dA24_diota*h24;
     dh[PSI]   = c->dA21_dpsi*h21  + c->dA22_dpsi*h22  + c->dA23_dpsi*h23  + c->dA24_dpsi*h24;
     dh[PHASE] = c->dA21_dphi0*h21 + c->dA22_dphi0*h22 + c->dA23_dphi0*h23 + c->dA24_dphi0*h24;
     dh[H0]    = c->dA21_dh0*h21   + c->dA22_dh0*h22   + c->dA23_dh0*h23   + c->dA24_dh0*h24;
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
 
 
// Fills the per-signal part of the context (parameters, amplitudes, phase constants)
static void setup_ctx(Deriv_ctx *c, const Fisher_settings *sett,
                      const Command_line_opts *opts, const Signal_params *sp)
{
     const double dt = sett->dt;
 
     c->alpha = sp->ra;
     c->delta = sp->dec;
     c->f     = sp->freq;
     c->fdot  = sp->fdot;
 
     amplitude_derivs(c, sp->h0, sp->iota, sp->psi, sp->phase);
 
     // phase constants, as in signal_gen()
     c->dt  = dt;
     c->K          = (double)sett->N * (opts->seg - sp->reffr);         // Nref
     c->freqo1     = M_PI*c->fdot*dt*dt;
     c->freqo0     = 2.*M_PI*c->f*dt + 2.*c->freqo1*c->K;
     c->phaseshift = c->freqo0*c->K - c->freqo1*c->K*c->K;
     c->wh         = 2.*M_PI*dt*sett->fpo;     // heterodyne, does not depend on the parameters
     c->omr = sett->omr;                // C_OMEGA_R * dt
}
 
// Fisher matrix of detector n:  F_jk = (dh_j | dh_k).
// If dA, dPsi are given, also  v_j = (Delta h | dh_j),  Delta h = dA h + dPsi dh/dPsi.
static void detector_sums(Deriv_ctx *c, int n, const Fisher_settings *sett,
                          const double *dA, const double *dPsi,
                          double F[NPARAMS][NPARAMS], double v[NPARAMS])
{
     const int N = sett->N;
 
     // rogcvir() coefficients and phi_r of this detector
     c->c1 = ifo[n].amod.c1;  c->c2 = ifo[n].amod.c2;  c->c3 = ifo[n].amod.c3;
     c->c4 = ifo[n].amod.c4;  c->c5 = ifo[n].amod.c5;  c->c6 = ifo[n].amod.c6;
     c->c7 = ifo[n].amod.c7;  c->c8 = ifo[n].amod.c8;  c->c9 = ifo[n].amod.c9;
     c->phir = ifo[n].sig.phir;
 
     for (int j = 0; j < NPARAMS; j++) {
          v[j] = 0.;
          for (int k = 0; k < NPARAMS; k++) F[j][k] = 0.;
     }
 
     for (int i = 0; i < N; i++) {
 
          if (!ifo[n].sig.noise[i]) continue;      // no signal in data gaps
 
          double h, dh[NPARAMS];
          signal_derivs(c, i, &ifo[n].sig.DetSSB[3*i], &h, dh);
 
          const double w = (i == 0 || i == N-1) ? 0.5 : 1.0;   // trapezoid
          for (int j = 0; j < NPARAMS; j++)
               for (int k = 0; k <= j; k++)
                    F[j][k] += w * dh[j] * dh[k];
 
          if (dA) {
               const double Dh = dA[i]*h + dPsi[i]*dh[PHASE];  // dh[PHASE] = dh/dPsi
               for (int j = 0; j < NPARAMS; j++)
                    v[j] += w * Dh * dh[j];
          }
     }
 
     const double norm = sett->dt / ifo[n].sig.var;
     for (int j = 0; j < NPARAMS; j++) {
          v[j] *= norm;
          for (int k = 0; k <= j; k++)
               F[k][j] = F[j][k] = norm * F[j][k];
     }
}
 
 
// C = F^-1 via the normalised matrix R = D^-1 F D^-1. Returns 0 on success.
static int covariance(const double F[NPARAMS][NPARAMS], double C[NPARAMS][NPARAMS])
{
     double d[NPARAMS], R[NPARAMS][NPARAMS], Ri[NPARAMS][NPARAMS];
 
     for (int j = 0; j < NPARAMS; j++) {
          if (!(F[j][j] > 0.)) return 1;
          d[j] = sqrt(F[j][j]);
     }
     for (int j = 0; j < NPARAMS; j++)
          for (int k = 0; k < NPARAMS; k++)
               R[j][k] = F[j][k] / (d[j]*d[k]);
 
     if (invm(&R[0][0], NPARAMS, &Ri[0][0])) return 1;
 
     for (int j = 0; j < NPARAMS; j++)
          for (int k = 0; k < NPARAMS; k++)
               C[j][k] = 0.5*(Ri[j][k] + Ri[k][j]) / (d[j]*d[k]);
     return 0;
}
 
 
static const char *pname[NPARAMS] = {"freq", "fdot", "ra", "dec", "iota", "psi", "phase", "h0"};
 
 
void fisher_matrix( Fisher_settings *sett,
                    Command_line_opts *opts,
                    Aux_arrays *aux_arr,
                    Signal_params *sig_params)
{
     // h0 from the SNR, if requested
     if (!strcmp(sig_params->amporsnr, "snr")) {
          signal_gen(sett, opts, aux_arr, sig_params, ifo);
          strcpy(sig_params->amporsnr, "amp");
     }
 
     Deriv_ctx c;
     setup_ctx(&c, sett, opts, sig_params);
 
     // network Fisher matrix = sum over detectors
     double F[NPARAMS][NPARAMS] = {{0}}, C[NPARAMS][NPARAMS];
     for (int n = 0; n < sett->nifo; n++) {
          double Fd[NPARAMS][NPARAMS], v[NPARAMS];
          detector_sums(&c, n, sett, NULL, NULL, Fd, v);
          for (int j = 0; j < NPARAMS; j++)
               for (int k = 0; k < NPARAMS; k++)
                    F[j][k] += Fd[j][k];
     }
 
     save_matrix(&F[0][0], NPARAMS, "Fisher matrix", "fisher_matrix", opts);
 
     if (covariance(F, C)) {
          fprintf(stderr, "fisher_matrix: Fisher matrix is singular, "
                          "covariance matrix not computed\n");
          return;
     }
     save_matrix(&C[0][0], NPARAMS, "Inverse Fisher matrix (covariance)",
                 "inv_fisher_matrix", opts);
 
     printf("\n1-sigma errors:\n");
     for (int j = 0; j < NPARAMS; j++)
          printf("  %-6s %.6e\n", pname[j], sqrt(C[j][j]));
 
} // fisher_matrix
 
 
void systematics( Fisher_settings *sett,
                  Command_line_opts *opts,
                  Aux_arrays *aux_arr,
                  Signal_params *sig_params)
{
     const int N = sett->N;
 
     if (!strcmp(sig_params->amporsnr, "snr")) {
          signal_gen(sett, opts, aux_arr, sig_params, ifo);
          strcpy(sig_params->amporsnr, "amp");
     }
 
     // dA and dPsi: two columns, one row per sample
     double *dA = malloc(N * sizeof *dA), *dPsi = malloc(N * sizeof *dPsi);
     FILE *in = fopen(opts->sysfile, "r");
     if (!in) {
          perror(opts->sysfile);
          exit(EXIT_FAILURE);
     }
     for (int i = 0; i < N; i++)
          if (fscanf(in, "%lf %lf", &dA[i], &dPsi[i]) != 2) {
               fprintf(stderr, "%s: could not read row %d (need N = %d rows)\n",
                       opts->sysfile, i+1, N);
               exit(EXIT_FAILURE);
          }
     fclose(in);
 
     char filename[FNAME_LENGTH + 64];
     snprintf(filename, sizeof(filename), "%s/systematics_%03d_%04d%s.txt",
              opts->outdir, opts->seg, opts->band, opts->label);
     FILE *out = fopen(filename, "w");
     if (!out) {
          perror(filename);
          exit(EXIT_FAILURE);
     }
     fprintf(out, "# detector param Delta_theta sigma Delta_theta/sigma\n");
 
     Deriv_ctx c;
     setup_ctx(&c, sett, opts, sig_params);
 
     for (int n = 0; n < sett->nifo; n++) {
 
          // Fisher matrix and (Delta h | dh_j) of this detector
          double F[NPARAMS][NPARAMS], C[NPARAMS][NPARAMS], v[NPARAMS];
          detector_sums(&c, n, sett, dA, dPsi, F, v);
 
          printf("\nSystematic errors, %.2s:\n", ifo[n].name);
          if (covariance(F, C)) {
               fprintf(stderr, "systematics: Fisher matrix of %.2s is singular\n", ifo[n].name);
               continue;
          }
 
          // Delta theta = F^-1 v
          printf("  %-6s %-14s %-14s %s\n", "param", "Delta_theta", "sigma", "Delta_theta/sigma");
          for (int j = 0; j < NPARAMS; j++) {
               double dth = 0.;
               for (int k = 0; k < NPARAMS; k++)
                    dth += C[j][k] * v[k];
               const double s = sqrt(C[j][j]);
               printf("  %-6s % .6e  %.6e  % .3e\n", pname[j], dth, s, dth/s);
               fprintf(out, "%.2s %s %.16e %.16e %.16e\n", ifo[n].name, pname[j], dth, s, dth/s);
          }
     }
 
     fclose(out);
     free(dA);
     free(dPsi);
 
} // systematics
