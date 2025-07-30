#ifndef __AUXI_H__
#define __AUXI_H__

#include <complex.h>

#define sqr(x) ((x)*(x))
#define TOSTRA(x) #x
#define TOSTR(x) TOSTRA(x)

#define TINY 1.0e-20
#define NINTERP 3  /* degree of the interpolation polynomial - do not change!!! */
#define NAVFSTAT 4096
#define round(x) floor((x)+0.5)

// Define COMP_FLOAT to switch to single precision of triggers and of fft
#ifdef COMP_FLOAT
    #define FLOAT_TYPE float
    #define FFTW_PRE(NAME) fftwf ## NAME
    #define NORM(x) ( crealf(x)*crealf(x) + cimagf(x)*cimagf(x) )
#else
    #define FLOAT_TYPE double
    #define FFTW_PRE(NAME) fftw ## NAME
    #define NORM(x) ( creal(x)*creal(x) + cimag(x)*cimag(x) )
#endif

void lin2ast(double be1, double be2, int pm, double sepsm, double cepsm,
             double *sinal, double *cosal, double *sindel, double *cosdel);

int ast2lin(double alfa, double delta, double epsm, double *be);
int dimless_to_code(double *M, double *pos_dimless, float *pos_code);

void spline(complex double *, int, complex double *);
complex double splint (complex double *, complex double *, int, double);
void splintpad (complex double *, double *, int, int, complex double*);
void linterp (complex double *, double *, int, int, complex double*);
void triginterp (complex double *ya, complex double *yb, double *shftf,
		       int N, int nfft, complex double *outa, complex double *outb);
double var (float *, int);

void gridr (double *, float *, float *, float *, double, double);
double FStat (FLOAT_TYPE *, int, int, int);

int ludcmp (double *, int, int *, double *);
int lubksb (double *, int, int *, double *);

// gridopt
int invm (const double *, int, double *);
double det (const double *, int);

// for qsorting the lines
int compared2c (const void *, const void *);

// signal handler
static void sig_handler(int signo);

#endif
