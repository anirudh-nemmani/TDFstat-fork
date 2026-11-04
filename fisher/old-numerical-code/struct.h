#ifndef __STRUCT_H__
#define __STRUCT_H__

#include <complex.h>
#include <hdf5.h>

#define MAX_DETECTORS 8        // Maximum number of detectors in network
#define DETNAME_LENGTH 2       // Detector name length (H1, L1, V1...)
#define FNAME_LENGTH 2048      // Maximum length of a filename

/* Command line options for Fisher */
typedef struct _comm_line_opts {

     int seg, band;
     double fpo_val, overlap;
     const char *indir, *outdir, *usedet, *label, *addsig, *mods;
     char state_file[FNAME_LENGTH];

} Command_line_opts;


/* Auxiluary arrays */
typedef struct _aux_arrays {

     double *sinmodf, *cosmodf; // Earth position
     double *t2;                // time^2

} Aux_arrays;

/* Fisher settings */
typedef struct _fisher_settings {

     double fpo,    // Band frequency
            dt,     // Sampling time
            B,      // Bandwidth
            oms,    // Dimensionless angular frequency (fpo)
            omr,    // C_OMEGA_R * dt (dimensionless Earth's angular frequency)
            sepsm,	// sin(epsm)
            cepsm;	// cos(epsm)

     int nod,          // number of days of observation
         N,            // number of data points
         s,            // number of spindowns
         nd,           // degrees of freedom
         nifo;         // number of detectors

} Fisher_settings;

/* Signal parameters */
typedef struct _signal_params {
    int reffr;                   // Reference frame number
    double h0, snr, freq, fdot;  // Intrinsic signal parameters
    double iota, psi, phase;     // Amplitude modulation parameters
    double ra, dec;              // sky position
    char amporsnr[4];                 // String of 3 characters + null terminator
} Signal_params;

/* input signal arrays */
typedef struct _signals {

     float *noise;         // Noise time series
     float *signal;        // Signal time series
     double *DetSSB;       // Ephemeris of the detector
     double *aa, *bb;      // Amplitude modulation functions
     double *shftf, *shft; // Resampling and time-shifting

     double epsm,
            phir,
            sepsm,	  // sin(epsm)
            cepsm,	  // cos(epsm)
            sphir,	  // sin(phi_r)
            cphir,	  // cos(phi_r)
            crf0,    // number of 0s as: N/(N-Nzeros)
            var; 	  // variance of the noise

     int Nzeros;
     double complex *xDatma, *xDatmb;

} Signals;

/* Amplitude modulation function coefficients */
typedef struct _ampl_mod_coeff {

	double c1, c2, c3, c4, c5, c6, c7, c8, c9;

} Ampl_mod_coeff;

/* Detector and its data related settings */
typedef struct _detector {

     char xdatname[FNAME_LENGTH];
     char name[DETNAME_LENGTH];

     double ephi, 		// Geographical latitude phi in radians
            elam, 		// Geographical longitude in radians
            eheight,     // Height h above the Earth ellipsoid in meters
            egam; 		// Orientation of the detector gamma

     Ampl_mod_coeff amod;
     Signals sig;

     double start_time;  // Start time of the detector

} Detector_settings;

/* Global array of detectors (network) */
extern Detector_settings ifo[MAX_DETECTORS];

typedef enum {
    FREQ,
    FDOT,
    RA,
    DEC,
    IOTA,
    PSI,
    PHASE
} ParamIndex;

#endif
