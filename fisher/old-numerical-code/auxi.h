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

void spline(complex double *, int, complex double *);
complex double splint (complex double *, complex double *, int, double);
void splintpad (complex double *, double *, int, int, complex double*);
void linterp (complex double *, double *, int, int, complex double*);
double var (float *, int);

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
