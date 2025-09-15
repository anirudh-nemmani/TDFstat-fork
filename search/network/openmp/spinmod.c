#define _GNU_SOURCE
#include <math.h>
#include <stdio.h>
#include <sys/types.h>
#include <complex.h>
#include <omp.h>

#include "auxi.h"
#include "settings.h"

#if defined(SLEEF)
#include <sleef.h>
//##define SLEEF_ENABLE_OMP_SIMD
#endif


void spindown_modulation(const int nifo, const int N, const FLOAT_TYPE het1,
                         const FLOAT_TYPE spnd, const FLOAT_TYPE _tmp1[][N],
                         FFTW_PRE(_complex) *fxa, FFTW_PRE(_complex) *fxb)
{
     int i, n;

     for (n=0; n<nifo; ++n) {

#if defined(SLEEF)
          // use simd sincos from the SLEEF library

#if defined(USE_AVX2)

#if defined(COMP_FLOAT)
          /******************************************************/
          /*         SINGLE PRECISION AVX2 version              */
          /******************************************************/
#define v_sincos Sleef_sincosf8_u35avx2
#define VLEN 8
          typedef __m256 vect;
          typedef Sleef___m256_2 vect2;

          int Nvec = (N/VLEN)*VLEN; // fully vectorizable range

#pragma omp parallel for schedule(static) default(shared) firstprivate(n)
          for (i=0; i<Nvec; i+=VLEN) {
               int j;
               vect2 v;
               vect a;
               a = _mm256_setr_ps(
                    het1*(i)   + spnd*_tmp1[n][i],
                    het1*(i+1) + spnd*_tmp1[n][i+1],
                    het1*(i+2) + spnd*_tmp1[n][i+2],
                    het1*(i+3) + spnd*_tmp1[n][i+3],
                    het1*(i+4) + spnd*_tmp1[n][i+4],
                    het1*(i+5) + spnd*_tmp1[n][i+5],
                    het1*(i+6) + spnd*_tmp1[n][i+6],
                    het1*(i+7) + spnd*_tmp1[n][i+7]
               );
               v = v_sincos(a);
               // exph = _c[j] - I*_p[j];
               // fxa[i+j] = ifo[0].sig.xDatma[i+j]*exph;
               vect mare = _mm256_setr_ps(
                    creal(ifo[n].sig.xDatma[i]),
                    creal(ifo[n].sig.xDatma[i+1]),
                    creal(ifo[n].sig.xDatma[i+2]),
                    creal(ifo[n].sig.xDatma[i+3]),
                    creal(ifo[n].sig.xDatma[i+4]),
                    creal(ifo[n].sig.xDatma[i+5]),
                    creal(ifo[n].sig.xDatma[i+6]),
                    creal(ifo[n].sig.xDatma[i+7])
               );
               vect maim = _mm256_setr_ps(
                    cimag(ifo[n].sig.xDatma[i]),
                    cimag(ifo[n].sig.xDatma[i+1]),
                    cimag(ifo[n].sig.xDatma[i+2]),
                    cimag(ifo[n].sig.xDatma[i+3]),
                    cimag(ifo[n].sig.xDatma[i+4]),
                    cimag(ifo[n].sig.xDatma[i+5]),
                    cimag(ifo[n].sig.xDatma[i+6]),
                    cimag(ifo[n].sig.xDatma[i+7])
               );
               vect vec1, vre, vim;
               vec1 = _mm256_mul_ps(mare, v.y);
               vre = _mm256_fmadd_ps(maim, v.x, vec1);
               vec1 = _mm256_mul_ps(mare, v.x);
               vim = _mm256_fmsub_ps(maim, v.y, vec1);
               if (n==0)
                    for(j=0; j<VLEN; ++j)
                         fxa[i+j] = ((float *)&vre)[j] + I*((float *)&vim)[j];
               else
                    for(j=0; j<VLEN; ++j)
                         fxa[i+j] += ((float *)&vre)[j] + I*((float *)&vim)[j];

               mare = _mm256_setr_ps(
                    creal(ifo[n].sig.xDatmb[i]),
                    creal(ifo[n].sig.xDatmb[i+1]),
                    creal(ifo[n].sig.xDatmb[i+2]),
                    creal(ifo[n].sig.xDatmb[i+3]),
                    creal(ifo[n].sig.xDatmb[i+4]),
                    creal(ifo[n].sig.xDatmb[i+5]),
                    creal(ifo[n].sig.xDatmb[i+6]),
                    creal(ifo[n].sig.xDatmb[i+7])
               );
               maim = _mm256_setr_ps(
                    cimag(ifo[n].sig.xDatmb[i]),
                    cimag(ifo[n].sig.xDatmb[i+1]),
                    cimag(ifo[n].sig.xDatmb[i+2]),
                    cimag(ifo[n].sig.xDatmb[i+3]),
                    cimag(ifo[n].sig.xDatmb[i+4]),
                    cimag(ifo[n].sig.xDatmb[i+5]),
                    cimag(ifo[n].sig.xDatmb[i+6]),
                    cimag(ifo[n].sig.xDatmb[i+7])
               );
               vec1 = _mm256_mul_ps(mare, v.y);
               vre = _mm256_fmadd_ps(maim, v.x, vec1);
               vec1 = _mm256_mul_ps(mare, v.x);
               vim = _mm256_fmsub_ps(maim, v.y, vec1);
               if (n==0)
                    for(j=0; j<VLEN; ++j)
                         fxb[i+j] = ((float *)&vre)[j] + I*((float *)&vim)[j];
               else
                    for(j=0; j<VLEN; ++j)
                         fxb[i+j] += ((float *)&vre)[j] + I*((float *)&vim)[j];
          } // for(i)


#else // COMP_FLOAT
          /******************************************************/
          /*         DOUBLE PRECISION AVX2 version              */
          /******************************************************/
#define v_sincos Sleef_sincosd4_u35avx2
#define VLEN 4
          typedef __m256d vect;
          typedef Sleef___m256d_2 vect2;

          int Nvec = (N/VLEN)*VLEN; // fully vectorizable range

#pragma omp parallel for schedule(static) default(shared) firstprivate(n)
          for (i=0; i<Nvec; i+=VLEN) {
               int j;
               vect2 v;
               vect a;
               a = _mm256_setr_pd(
                    het1*(i)   + spnd*_tmp1[n][i],
                    het1*(i+1) + spnd*_tmp1[n][i+1],
                    het1*(i+2) + spnd*_tmp1[n][i+2],
                    het1*(i+3) + spnd*_tmp1[n][i+3]
               );
               v = v_sincos(a);
               // exph = _c[j] - I*_p[j];
               // fxa[i+j] = ifo[0].sig.xDatma[i+j]*exph;
               vect mare = _mm256_setr_pd(
                    creal(ifo[n].sig.xDatma[i]),
                    creal(ifo[n].sig.xDatma[i+1]),
                    creal(ifo[n].sig.xDatma[i+2]),
                    creal(ifo[n].sig.xDatma[i+3])
               );
               vect maim = _mm256_setr_pd(
                    cimag(ifo[n].sig.xDatma[i]),
                    cimag(ifo[n].sig.xDatma[i+1]),
                    cimag(ifo[n].sig.xDatma[i+2]),
                    cimag(ifo[n].sig.xDatma[i+3])
               );
               vect vec1, vre, vim;
               vec1 = _mm256_mul_pd(mare, v.y);
               vre = _mm256_fmadd_pd(maim, v.x, vec1);
               vec1 = _mm256_mul_pd(mare, v.x);
               vim = _mm256_fmsub_pd(maim, v.y, vec1);
               if (n==0)
                    for(j=0; j<VLEN; ++j)
                         fxa[i+j] = ((double *)&vre)[j] + I*((double *)&vim)[j];
               else
                    for(j=0; j<VLEN; ++j)
                         fxa[i+j] += ((double *)&vre)[j] + I*((double *)&vim)[j];

               mare = _mm256_setr_pd(
                    creal(ifo[n].sig.xDatmb[i]),
                    creal(ifo[n].sig.xDatmb[i+1]),
                    creal(ifo[n].sig.xDatmb[i+2]),
                    creal(ifo[n].sig.xDatmb[i+3])
               );
               maim = _mm256_setr_pd(
                    cimag(ifo[n].sig.xDatmb[i]),
                    cimag(ifo[n].sig.xDatmb[i+1]),
                    cimag(ifo[n].sig.xDatmb[i+2]),
                    cimag(ifo[n].sig.xDatmb[i+3])
               );
               vec1 = _mm256_mul_pd(mare, v.y);
               vre = _mm256_fmadd_pd(maim, v.x, vec1);
               vec1 = _mm256_mul_pd(mare, v.x);
               vim = _mm256_fmsub_pd(maim, v.y, vec1);
               if (n==0)
                    for(j=0; j<VLEN; ++j)
                         fxb[i+j] = ((double *)&vre)[j] + I*((double *)&vim)[j];
               else
                    for(j=0; j<VLEN; ++j)
                         fxb[i+j] += ((double *)&vre)[j] + I*((double *)&vim)[j];
          } // i

#endif //COMP_FLOAT

#elif defined(USE_AVX512)

#if defined(COMP_FLOAT)
          /******************************************************/
          /*       SINGLE PRECISION AVX512 version              */
          /******************************************************/
#define v_sincos Sleef_sincosf16_u35avx512f
#define VLEN 16
          typedef __m512 vect;
          typedef Sleef___m512_2 vect2;

          int Nvec = (N/VLEN)*VLEN; // fully vectorizable range

#pragma omp parallel for schedule(static) default(shared) firstprivate(n)
          for (i=0; i<Nvec; i+=VLEN) {
               int j;
               vect2 v;
               vect a;
               a = _mm512_setr_ps(
                    het1*(i)    + spnd*_tmp1[n][i],
                    het1*(i+1)  + spnd*_tmp1[n][i+1],
                    het1*(i+2)  + spnd*_tmp1[n][i+2],
                    het1*(i+3)  + spnd*_tmp1[n][i+3],
                    het1*(i+4)  + spnd*_tmp1[n][i+4],
                    het1*(i+5)  + spnd*_tmp1[n][i+5],
                    het1*(i+6)  + spnd*_tmp1[n][i+6],
                    het1*(i+7)  + spnd*_tmp1[n][i+7],
                    het1*(i+8)  + spnd*_tmp1[n][i+8],
                    het1*(i+9)  + spnd*_tmp1[n][i+9],
                    het1*(i+10) + spnd*_tmp1[n][i+10],
                    het1*(i+11) + spnd*_tmp1[n][i+11],
                    het1*(i+12) + spnd*_tmp1[n][i+12],
                    het1*(i+13) + spnd*_tmp1[n][i+13],
                    het1*(i+14) + spnd*_tmp1[n][i+14],
                    het1*(i+15) + spnd*_tmp1[n][i+15]
               );
               v = v_sincos(a);
               // exph = _c[j] - I*_p[j];
               // fxa[i+j] = ifo[0].sig.xDatma[i+j]*exph;
               vect mare = _mm512_setr_ps(
                    creal(ifo[n].sig.xDatma[i]),
                    creal(ifo[n].sig.xDatma[i+1]),
                    creal(ifo[n].sig.xDatma[i+2]),
                    creal(ifo[n].sig.xDatma[i+3]),
                    creal(ifo[n].sig.xDatma[i+4]),
                    creal(ifo[n].sig.xDatma[i+5]),
                    creal(ifo[n].sig.xDatma[i+6]),
                    creal(ifo[n].sig.xDatma[i+7]),
                    creal(ifo[n].sig.xDatma[i+8]),
                    creal(ifo[n].sig.xDatma[i+9]),
                    creal(ifo[n].sig.xDatma[i+10]),
                    creal(ifo[n].sig.xDatma[i+11]),
                    creal(ifo[n].sig.xDatma[i+12]),
                    creal(ifo[n].sig.xDatma[i+13]),
                    creal(ifo[n].sig.xDatma[i+14]),
                    creal(ifo[n].sig.xDatma[i+15])
               );
               vect maim = _mm512_setr_ps(
                    cimag(ifo[n].sig.xDatma[i]),
                    cimag(ifo[n].sig.xDatma[i+1]),
                    cimag(ifo[n].sig.xDatma[i+2]),
                    cimag(ifo[n].sig.xDatma[i+3]),
                    cimag(ifo[n].sig.xDatma[i+4]),
                    cimag(ifo[n].sig.xDatma[i+5]),
                    cimag(ifo[n].sig.xDatma[i+6]),
                    cimag(ifo[n].sig.xDatma[i+7]),
                    cimag(ifo[n].sig.xDatma[i+8]),
                    cimag(ifo[n].sig.xDatma[i+9]),
                    cimag(ifo[n].sig.xDatma[i+10]),
                    cimag(ifo[n].sig.xDatma[i+11]),
                    cimag(ifo[n].sig.xDatma[i+12]),
                    cimag(ifo[n].sig.xDatma[i+13]),
                    cimag(ifo[n].sig.xDatma[i+14]),
                    cimag(ifo[n].sig.xDatma[i+15])
               );
               vect vec1, vre, vim;
               vec1 = _mm512_mul_ps(mare, v.y);
               vre = _mm512_fmadd_ps(maim, v.x, vec1);
               vec1 = _mm512_mul_ps(mare, v.x);
               vim = _mm512_fmsub_ps(maim, v.y, vec1);
               if (n==0)
                    for(j=0; j<VLEN; ++j)
                         fxa[i+j] = ((float *)&vre)[j] + I*((float *)&vim)[j];
               else
                    for(j=0; j<VLEN; ++j)
                         fxa[i+j] += ((float *)&vre)[j] + I*((float *)&vim)[j];

               mare = _mm512_setr_ps(
                    creal(ifo[n].sig.xDatmb[i]),
                    creal(ifo[n].sig.xDatmb[i+1]),
                    creal(ifo[n].sig.xDatmb[i+2]),
                    creal(ifo[n].sig.xDatmb[i+3]),
                    creal(ifo[n].sig.xDatmb[i+4]),
                    creal(ifo[n].sig.xDatmb[i+5]),
                    creal(ifo[n].sig.xDatmb[i+6]),
                    creal(ifo[n].sig.xDatmb[i+7]),
                    creal(ifo[n].sig.xDatmb[i+8]),
                    creal(ifo[n].sig.xDatmb[i+9]),
                    creal(ifo[n].sig.xDatmb[i+10]),
                    creal(ifo[n].sig.xDatmb[i+11]),
                    creal(ifo[n].sig.xDatmb[i+12]),
                    creal(ifo[n].sig.xDatmb[i+13]),
                    creal(ifo[n].sig.xDatmb[i+14]),
                    creal(ifo[n].sig.xDatmb[i+15])
               );
               maim = _mm512_setr_ps(
                    cimag(ifo[n].sig.xDatmb[i]),
                    cimag(ifo[n].sig.xDatmb[i+1]),
                    cimag(ifo[n].sig.xDatmb[i+2]),
                    cimag(ifo[n].sig.xDatmb[i+3]),
                    cimag(ifo[n].sig.xDatmb[i+4]),
                    cimag(ifo[n].sig.xDatmb[i+5]),
                    cimag(ifo[n].sig.xDatmb[i+6]),
                    cimag(ifo[n].sig.xDatmb[i+7]),
                    cimag(ifo[n].sig.xDatmb[i+8]),
                    cimag(ifo[n].sig.xDatmb[i+9]),
                    cimag(ifo[n].sig.xDatmb[i+10]),
                    cimag(ifo[n].sig.xDatmb[i+11]),
                    cimag(ifo[n].sig.xDatmb[i+12]),
                    cimag(ifo[n].sig.xDatmb[i+13]),
                    cimag(ifo[n].sig.xDatmb[i+14]),
                    cimag(ifo[n].sig.xDatmb[i+15])
               );
               vec1 = _mm512_mul_ps(mare, v.y);
               vre = _mm512_fmadd_ps(maim, v.x, vec1);
               vec1 = _mm512_mul_ps(mare, v.x);
               vim = _mm512_fmsub_ps(maim, v.y, vec1);
               if (n==0)
                    for(j=0; j<VLEN; ++j)
                         fxb[i+j] = ((float *)&vre)[j] + I*((float *)&vim)[j];
               else
                    for(j=0; j<VLEN; ++j)
                         fxb[i+j] += ((float *)&vre)[j] + I*((float *)&vim)[j];
          } // i


#else // COMP_FLOAT
          /******************************************************/
          /*       DOUBLE PRECISION AVX512 version              */
          /******************************************************/
#define v_sincos Sleef_sincosd8_u35avx512f
#define VLEN 8
          typedef __m512d vect;
          typedef Sleef___m512d_2 vect2;

          int Nvec = (N/VLEN)*VLEN; // fully vectorizable range

#pragma omp parallel for schedule(static) default(shared) firstprivate(n)
          for (i=0; i<Nvec; i+=VLEN) {
               int j;
               vect2 v;
               vect a;
               a = _mm512_setr_pd(
                    het1*(i)   + spnd*_tmp1[n][i],
                    het1*(i+1) + spnd*_tmp1[n][i+1],
                    het1*(i+2) + spnd*_tmp1[n][i+2],
                    het1*(i+3) + spnd*_tmp1[n][i+3],
                    het1*(i+4) + spnd*_tmp1[n][i+4],
                    het1*(i+5) + spnd*_tmp1[n][i+5],
                    het1*(i+6) + spnd*_tmp1[n][i+6],
                    het1*(i+7) + spnd*_tmp1[n][i+7]
               );
               v = v_sincos(a);
               // exph = _c[j] - I*_p[j];
               // fxa[i+j] = ifo[0].sig.xDatma[i+j]*exph;
               vect mare = _mm512_setr_pd(
                    creal(ifo[n].sig.xDatma[i]),
                    creal(ifo[n].sig.xDatma[i+1]),
                    creal(ifo[n].sig.xDatma[i+2]),
                    creal(ifo[n].sig.xDatma[i+3]),
                    creal(ifo[n].sig.xDatma[i+4]),
                    creal(ifo[n].sig.xDatma[i+5]),
                    creal(ifo[n].sig.xDatma[i+6]),
                    creal(ifo[n].sig.xDatma[i+7])
               );
               vect maim = _mm512_setr_pd(
                    cimag(ifo[n].sig.xDatma[i]),
                    cimag(ifo[n].sig.xDatma[i+1]),
                    cimag(ifo[n].sig.xDatma[i+2]),
                    cimag(ifo[n].sig.xDatma[i+3]),
                    cimag(ifo[n].sig.xDatma[i+4]),
                    cimag(ifo[n].sig.xDatma[i+5]),
                    cimag(ifo[n].sig.xDatma[i+6]),
                    cimag(ifo[n].sig.xDatma[i+7])
               );
               vect vec1, vre, vim;
               vec1 = _mm512_mul_pd(mare, v.y);
               vre = _mm512_fmadd_pd(maim, v.x, vec1);
               vec1 = _mm512_mul_pd(mare, v.x);
               vim = _mm512_fmsub_pd(maim, v.y, vec1);
               if (n==0)
                    for(j=0; j<VLEN; ++j)
                         fxa[i+j] = ((double *)&vre)[j] + I*((double *)&vim)[j];
               else
                    for(j=0; j<VLEN; ++j)
                         fxa[i+j] += ((double *)&vre)[j] + I*((double *)&vim)[j];

               mare = _mm512_setr_pd(
                    creal(ifo[n].sig.xDatmb[i]),
                    creal(ifo[n].sig.xDatmb[i+1]),
                    creal(ifo[n].sig.xDatmb[i+2]),
                    creal(ifo[n].sig.xDatmb[i+3]),
                    creal(ifo[n].sig.xDatmb[i+4]),
                    creal(ifo[n].sig.xDatmb[i+5]),
                    creal(ifo[n].sig.xDatmb[i+6]),
                    creal(ifo[n].sig.xDatmb[i+7])
               );
               maim = _mm512_setr_pd(
                    cimag(ifo[n].sig.xDatmb[i]),
                    cimag(ifo[n].sig.xDatmb[i+1]),
                    cimag(ifo[n].sig.xDatmb[i+2]),
                    cimag(ifo[n].sig.xDatmb[i+3]),
                    cimag(ifo[n].sig.xDatmb[i+4]),
                    cimag(ifo[n].sig.xDatmb[i+5]),
                    cimag(ifo[n].sig.xDatmb[i+6]),
                    cimag(ifo[n].sig.xDatmb[i+7])
               );
               vec1 = _mm512_mul_pd(mare, v.y);
               vre = _mm512_fmadd_pd(maim, v.x, vec1);
               vec1 = _mm512_mul_pd(mare, v.x);
               vim = _mm512_fmsub_pd(maim, v.y, vec1);
               if (n==0)
                    for(j=0; j<VLEN; ++j)
                         fxb[i+j] = ((double *)&vre)[j] + I*((double *)&vim)[j];
               else
                    for(j=0; j<VLEN; ++j)
                         fxb[i+j] += ((double *)&vre)[j] + I*((double *)&vim)[j];
          } // i

#endif // COMP_FLOAT

#elif defined(USE_AVX2_SINCOS)
          // use only SIMD sincos from sleef
          // only double version

#define v_sincos Sleef_sincosd4_u35
#define VLEN 4
          typedef __m256d vect;
          typedef Sleef___m256d_2 vect2;
          static vect vload_vd_p(const double *ptr) { return _mm256_load_pd(ptr); }
          static void vstore_v_p_vd(double *ptr, vect v) { _mm256_store_pd(ptr, v); }

#pragma omp parallel for schedule(static) default(shared) firstprivate(n)
          for (i=0; i<Nvec; i+=VLEN) {
               int j;
               double _p[VLEN] __attribute__((aligned(32))),
                      _c[VLEN] __attribute__((aligned(32)));
               complex double exph;
               vect2 v;
               vect a;

               for(j=0; j<VLEN; ++j)
                    _p[j] =  het1*(i+j) + spnd*_tmp1[n][i+j];
               a = vload_vd_p(_p);
               v = v_sincos(a);
               vstore_v_p_vd(_p, v.x); // reuse _p for sin
               vstore_v_p_vd(_c, v.y);

               if (n==0)
                    for(j=0; j<VLEN; ++j){
                         exph = _c[j] - I*_p[j];
                         fxa[i+j] = ifo[n].sig.xDatma[i+j]*exph;
                         fxb[i+j] = ifo[n].sig.xDatmb[i+j]*exph;
                    }
               else
                    for(j=0; j<VLEN; ++j){
                         exph = _c[j] - I*_p[j];
                         fxa[i+j] += ifo[n].sig.xDatma[i+j]*exph;
                         fxb[i+j] += ifo[n].sig.xDatmb[i+j]*exph;
                    }
          } // i
#endif // USE_AVX2_SINCOS

          /*****************************************************************/
          /*  all vectorized versions, if the last vector is incoplete     */
          /*****************************************************************/
          for(i=Nvec; i<N; ++i) {
               double phase, sp, cp;
               complex double exph;
               phase = het1*i + spnd*_tmp1[n][i];
               sincos(phase, &sp, &cp);
               exph = cp - I*sp;
               if (n==0){
                    fxa[i] = ifo[n].sig.xDatma[i]*exph;
                    fxb[i] = ifo[n].sig.xDatmb[i]*exph;
               } else {
                    fxa[i] += ifo[n].sig.xDatma[i]*exph;
                    fxb[i] += ifo[n].sig.xDatmb[i]*exph;
               }
          } // i

#elif defined(GNUSINCOS)  // !SLEEF
#if defined(COMP_FLOAT)
#define sincos sincosf
#else
#define sincos sincos
#endif
#pragma omp parallel for schedule(static) default(shared) firstprivate(n)
          for(i=N-1; i>-1; --i) {
               FLOAT_TYPE phase, cp, sp;
               complex FLOAT_TYPE exph;

               phase = het1*i + spnd*_tmp1[n][i];
               // full phase formula test
               // phase += spnd * ifo[n].sig.shft[i]*ifo[n].sig.shft[i];
               sincos(phase, &sp, &cp);
               exph = cp - I*sp;
               if (n==0){
                    fxa[i] = ifo[n].sig.xDatma[i]*exph;
                    fxb[i] = ifo[n].sig.xDatmb[i]*exph;
               } else {
                    fxa[i] += ifo[n].sig.xDatma[i]*exph;
                    fxb[i] += ifo[n].sig.xDatmb[i]*exph;
               }

          }
#else
     // generic version for systems without sincos function
#pragma omp parallel for schedule(static) default(shared) firstprivate(n)
          for(i=N-1; i!=-1; --i) {
               double phase, cp, sp;
               complex double exph;

               phase = het1*i + spnd*_tmp1[n][i];
               cp = cos(phase);
               sp = sin(phase);
               exph = cp - I*sp;
               if (n==0){
                    fxa[i] = ifo[n].sig.xDatma[i]*exph;
                    fxb[i] = ifo[n].sig.xDatmb[i]*exph;
               } else {
                    fxa[i] += ifo[n].sig.xDatma[i]*exph;
                    fxb[i] += ifo[n].sig.xDatmb[i]*exph;
               }
          }
#endif

     } // n

} // spindown_modulation
