#ifndef __SETTINGS_H__
#define __SETTINGS_H__

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <complex.h>
#include <float.h>

#include "auxi.h"

#include "gsl/gsl_vector.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#include <omp.h>


inline void spline(complex double *y, int n, complex double *y2)
{
     int i, k;
     complex double invp, qn, un;
     static complex double *u = NULL;

     if (!u) u = (complex double *)malloc((n-1)*sizeof(complex double));
     y2[0] = u[0] = 0.;

     for (i=1; i<n-1; ++i) {
          //p = .5*y2[i-1]+2.;
          //y2[i] = -.5/p;
          //u[i] = y[i+1]-2.*y[i]+y[i-1];
          //u[i] = (3.*u[i]-.5*u[i-1])/p;
          invp = 2./(y2[i-1]+4.);
          y2[i] = -.5*invp;
          u[i] = y[i-1]-2.*y[i]+y[i+1];
          u[i] = (-.5*u[i-1]+3.*u[i])*invp;
     }
     qn = un = 0.;
     y2[n-1] = (un-qn*u[n-2])/(qn*y2[n-2]+1.);
     for (k=n-2; k>=0; --k)
          y2[k] = y2[k]*y2[k+1]+u[k];

} /* spline() */


inline complex double splint (complex double *ya, complex double *y2a, int n, double x)
{
     int klo, khi;
     double b, a;

     if (x<0 || x>n-1)
          return 0.;
     klo = floor (x);
     khi = klo+1;
     a = khi - x;
     b = x - klo;
     return a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])/6.0;

} /* splint() */


void splintpad (complex double *ya, double *shftf, int N, int interpftpad,
                complex double *out)
{
     /* Cubic spline with "natural" boundary conditions.
     Input:
     ya[i] - value of the function being interpolated in x_i = i,
     for i = 0 .. (interpftpad*N-1)	(changed on exit);
     Interpolating spline will be calculated at the points
     interpftpad*(i-shftf[i]), for i = 0 .. (N-1);
     N - number of output data points.
     Output:
     out[i] - value of the interpolating function
     at interpftpad*(i-shftf[i]).
     */
     complex double *y2;
     double x;
     int i;

     y2 = (complex double *) malloc (interpftpad*N*sizeof (complex double)); //vector twice-size of N
     spline (ya, interpftpad*N, y2);
#if defined(_OPENMP)
#pragma omp parallel default(shared) private(x)
#endif
     {
#if defined(_OPENMP)
#pragma omp for schedule(static)
#endif
          for (i=0; i<N; ++i) {
               x = interpftpad*(i-shftf[i]);
               out[i] = splint (ya, y2, interpftpad*N, x);
          } /* for i */
     }
     free (y2);
} /* splintpad */


// pci test
void linterp (complex double *ya, double *shftf, int N, int interpftpad, complex double *out)
{
     /* linear interpolation
     Input:
	ya[i] - value of the function being interpolated in x_i = i,
	for i = 0 .. (interpftpad*N-1)
	Interpolating spline will be calculated at the points interpftpad*(i-shftf[i]), for i = 0 .. (N-1);
	N - number of output data points.
	Output:
	out[i] - value of the interpolating function at interpftpad*(i-shftf[i]).
     */

     double x;
     int i;

     for (i=0; i<N; ++i) {
          x = interpftpad*(i-shftf[i]);
          int i1 = (int) x;
          int i2 = i1+1;
          double dx = x-i1;
          double dya_abs = cabs(ya[i2]) - cabs(ya[i1]);
          double dya_arg = carg(ya[i2]) - carg(ya[i1]);
          //out[i] = (ya[i2]-ya[i1])*dx;
          out[i] = dya_abs*dx*cexp(dya_arg*dx*I);
     }

} //linterp()


double var (float *x, int n)
{
     /* var(x, n) returns the variance (square of the standard deviation)
     of a given vector x of length n.
     */
     int i;
     double mean=0., variance=0.;

     for (i=0; i<n; i++)
          mean += x[i];
     mean /= n;
     for (i=0; i<n; i++)
          variance += sqr (x[i]-mean);
     variance /= (n-1);
     return variance;
} /* var() */


int ludcmp (double *a, int n, int *indx, double *d)
{
     /*	LU decomposition of a given real matrix a[0..n-1][0..n-1]
     Input:
     a		- an array containing elements of matrix a
                 (changed on exit)
     n		- number of rows and columns of a
     Output:
     indx - row permutation effected by the partial pivoting
     d	- +-1 depending on whether the number of rows
     interchanged was even or odd, respectively
     */

     int i, imax = -1, j, k;
     double big, dum, sum, temp;
     double *vv;

     vv = (double *) calloc (n, sizeof (double));
     *d = 1.0;
     for (i=0; i<n; i++) {
          big = 0.0;
          for (j=0; j<n; j++)
               if ((temp=fabs (a[n*i+j])) > big)
                    big = temp;
          if (big == 0.0)
               return 1;
          vv[i] = 1.0/big;
     }
     for (j=0; j<n; j++) {
          for (i=0; i<j; i++) {
               sum = a[n*i+j];
               for (k=0; k<i; k++)
                    sum -= a[n*i+k]*a[n*k+j];
               a[n*i+j] = sum;
          }
          big = 0.0;
          for (i=j; i<n; i++) {
               sum = a[n*i+j];
               for (k=0; k<j; k++)
                    sum -= a[n*i+k]*a[n*k+j];
               a[n*i+j] = sum;
               if ((dum = vv[i]*fabs (sum)) >= big) {
                    big = dum;
                    imax = i;
               }
          }
          if (j != imax) {
               for (k=0; k<n; k++) {
                    dum = a[n*imax+k];
                    a[n*imax+k] = a[n*j+k];
                    a[n*j+k] = dum;
               }
               *d = -(*d);
               vv[imax] = vv[j];
          }
          indx[j] = imax;
          if (a[n*j+j] == 0.0)
               a[n*j+j] = TINY;
          if (j != n) {
               dum = 1.0/(a[n*j+j]);
               for (i=j+1; i<n; i++)
                    a[n*i+j] *= dum;
          }
     }
     free (vv);
     return 0;
} /* ludcmp() */


int lubksb (double *a, int n, int *indx, double *b)
{
     /* Solves the set of n linear equations A X=B.
     Input:
     a[0..n-1][0..n-1] - LU decomposition af a matrix A,
     determined by ludcmp()
     n				- number of rows and columns of a
     indx[0..n-1]		- permutation vector returned by ludcmp
     b[0..n-1]			 - right-hand side vector B
     (changed on exit)
     Output:
     b[0..n-1]			- solution vector X
     */

     int i, ii=-1, ip, j;
     double sum;

     for (i=0; i<n; i++) {
          ip = indx[i];
          sum = b[ip];
          b[ip] = b[i];
          if (ii>=0)
               for (j=ii; j<=i-1; j++)
                    sum -= a[n*i+j]*b[j];
          else if (sum)
               ii = i;
          b[i] = sum;
     }
     for (i=n-1; i>=0; i--) {
          sum = b[i];
          for (j=i+1; j<n; j++)
               sum -= a[n*i+j]*b[j];
          b[i] = sum/a[n*i+i];
     }
     return 0;
} /* lubksb() */


int invm (const double *a, int N, double *y)
{
     /* Inverse of a real matrix a[0..N-1][0..N-1].
     Input:
     a[0..N-1][0..N-1] - given matrix (saved on exit)
     N	      - number of rows and columns of a
     Output:
     y[0..N-1][0..N-1] - inverse of a
     */

     double d, *col, *al;
     int i, j, *indx;

     al = (double *) calloc (sqr(N), sizeof (double));
     indx = (int *) calloc (N, sizeof (int));
     col = (double *) calloc (N, sizeof (double));
     for (i=0; i<sqr(N); i++)
          al[i] = a[i];
     if (ludcmp (al, N, indx, &d))
          return 1;
     for (j=0; j<N; j++) {
          for (i=0; i<N; i++)
               col[i] = 0.0;
          col[j] = 1.0;
          lubksb (al, N, indx, col);
          for (i=0; i<N; i++)
               y[N*i+j] = col[i];
     }
     free (col);
     free (indx);
     free (al);
     return 0;
} /* invm() */


double det (const double *a, int N)
{
     /* determinant of a real matrix a[0..N-1][0..N-1] */
     double d, *al;
     int j, *indx;

     al = (double *) calloc (sqr(N), sizeof (double));
     indx = (int *) calloc (N, sizeof (int));
     for (j=0; j<sqr(N); j++)
          al[j] = a[j];
     ludcmp (al, N, indx, &d);
     for (j=0; j<N; j++)
          d *= al[N*j+j];
     free (indx);
     free (al);
     return d;
} /* det() */


int compared2c(const void *a, const void *b)
{
     double* da = (double*)a;
     double* db = (double*)b;

     int diff1 = (da[0] > db[0]) - (da[0] < db[0]);
     if (diff1 != 0) return diff1;
     return (da[1] > db[1]) - (da[1] < db[1]);
}

#endif
