#ifndef __COMMON_H__
#define __COMMON_H__

// auxi.h first: struct.h uses the FFTW_PRE macro defined there
#include "auxi.h"
#include "struct.h"

#define RAD_TO_DEG (180/M_PI) // = 180/pi

//constants
#define C_SPEED_OF_LIGHT 299792.458 // in km/s
#define C_AU 1.49597870691e8	// Astronomical unit, km

#define C_EPSMA (84381.448/3600./RAD_TO_DEG)
//#define C_EPSMA 0.409092804222328965124688693322241306304931640625
// Average obliquity of the ecliptic: 23.439

//Earth ellipsoid
#define C_ELLIPSOID_A 6378.140
#define C_ELLIPSOID_F 298.257
#define C_ELLIPSOID_B (C_ELLIPSOID_A * ( 1. - 1. / C_ELLIPSOID_F ) )
//6356.755288157528

#define C_OMEGA_R 7.2921151467064e-5
#define C_SIDDAY (2.*M_PI/C_OMEGA_R)
// 86164.09890369719 // 2.*M_PI/Omegar      // Sideral day
#define C_TAIDAY  86400.				            // TAI day

#define C_YEARSEC (365.25*C_TAIDAY)
//31557600.0          // year in seconds = 365.25 * 86400

double fpo(int band, double overlap, double dt);
void get_fdot_range(double fpo, double dt, double *fdot_min, double *fdot_max);

// Per-detector input file readers; det is an index into the ifo[] array.
// Called for every detector by init_arrays() (init.c), and individually by
// tools that need only part of the input data (e.g. vlinegen).
void read_xdat( Search_settings *sett, Command_line_opts *opts, int det );
void read_detssb( Search_settings *sett, Command_line_opts *opts, int det );
void read_start_time( Search_settings *sett, Command_line_opts *opts, int det );

#endif
