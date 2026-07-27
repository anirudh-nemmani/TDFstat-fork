// general purpose functions for the project

#include "common.h"

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
