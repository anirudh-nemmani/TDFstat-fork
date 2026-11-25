#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
//#include <getopt.h>
#include <gsl/gsl_linalg.h>
#include <time.h>
#include <dirent.h>
#include <signal.h>

#include "auxi.h"
#include "struct.h"
#include "common.h"
#include "settings.h"
#include "jobcore.h"
#include "init.h"

#ifndef CODEVER
#define CODEVER unknown
#endif

Detector_settings ifo[MAX_DETECTORS];
volatile sig_atomic_t save_state = 0;


int main (int argc, char* argv[])
{

    Command_line_opts opts;
    Search_settings sett;
    Search_range s_range;
    Aux_arrays aux_arr;

    // WIP INJ
    Signal_params sgnl_params;  // to be removed
    Injection_params *inj=NULL;
    int inj_idx=0; // default for inj_mode=multi

    int i;

    printf("git commit: %s\n", CODEVER);
    if (signal(SIGUSR1, sig_handler) != SIG_ERR &&
        signal(SIGTERM, sig_handler) != SIG_ERR ){
            printf("State saved on SIGTERM or SIGUSR1\n");
    }
    // Command line options
    read_ini_file(&sett, &opts, argc, argv);

    // Output data handling
    struct stat buffer;

    if (stat(opts.outdir, &buffer) == -1) {
        if (errno == ENOENT) {
            // Output directory apparently does not exist, try to create one
            if(mkdir(opts.outdir, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH	| S_IXOTH) == -1) {
                perror (opts.outdir);
                return 1;
            }
        } else { // can't access output directory
            perror (opts.outdir);
            return 1;
        }
    }

    // Detector network settings
    detectors_settings(&sett, &opts);

    // Grid data
    read_grid(&sett, &opts);

    // Search settings
    search_settings(&sett);

    // Array initialization and reading the ephemerids
    init_arrays(&sett, &opts, &aux_arr);

    // Narrowing-down the band (excluding the edges
    // according to the opts.narrowdown parameter)
    // adds two lines
    if(opts.narrowdown < 0.5*M_PI) narrow_down_band(&sett, &opts);

    // Handle veto lines and opts.veto_flag
    char band_vl_file[FNAME_LENGTH];
    sprintf(band_vl_file, "%s/vlines_%04d%s.dat", opts.outdir, opts.band, opts.label);
    if (access(band_vl_file, R_OK) == 0) {
        // read lines from band_vl_file
        read_band_vlines(&sett, &opts, band_vl_file);
    } else {
        // extract inband veto lines from global veto file
        // and write them to band_vl_file
        extract_band_vlines(&sett, &opts, band_vl_file);
    }

    // Amplitude modulation functions for each detector
    for(i=0; i<sett.nifo; i++)
        rogcvir(&ifo[i]);

    // Read injection file/dict to add injections (if multi mode)
    int ninj = 0;
    if (opts.inj_mode) {
        ninj = read_injection(&opts, &inj);
        if (!strcmp(opts.inj_mode, "multi"))
            add_injection(&sett, &opts, &aux_arr, inj, ninj, inj_idx);
    }
    // gcenter = injection: grid centered at the injection (only with inj_mode = serial, checked in read_ini_file)
    int gc_inj = strcmp(opts.gtype, "allsky") != 0 && opts.gcenter && !strcmp(opts.gcenter, "injection");

     // establish search range (s_range); gcenter = injection: centered at the first injection
     set_search_range(&sett, &opts, &s_range, gc_inj ? &inj[0] : NULL);

     // FFT plans
     FFTW_plans fftw_plans;
     FFTW_arrays fftw_arr;
     plan_fftw(&sett, &opts, &fftw_plans, &fftw_arr, &aux_arr);

     // Checkpointing
     int Fnum=0;	// candidate signal number
     read_checkpoints(&opts, &s_range, &Fnum);

    // main search job
    if (opts.inj_mode && !strcmp(opts.inj_mode, "serial")) {
        for (inj_idx=0; inj_idx<ninj; inj_idx++) {
            // gcenter = injection: grid centered at the current injection, search from its start
            if (gc_inj && inj_idx > 0) {
                set_search_range(&sett, &opts, &s_range, &inj[inj_idx]);
                // WIP: small steps to bypass read_checkpoints for injections
                s_range.pst = s_range.pmr[0];
                s_range.mst = s_range.mr[0];
                s_range.nst = s_range.nr[0];
                s_range.sst = s_range.spndr[0];
            }
            // add injections (if serial mode)
            add_injection(&sett, &opts, &aux_arr, inj, ninj, inj_idx);
            search(&sett, &opts, &s_range, &fftw_plans, &fftw_arr, &aux_arr, &Fnum, inj, ninj, inj_idx);
        }
    } else {
        // no injections or multi mode (just one search)
        search(&sett, &opts, &s_range, &fftw_plans, &fftw_arr, &aux_arr, &Fnum, inj, ninj, inj_idx);
    }

     // Cleanup & memory free
     cleanup(&sett, &opts, &s_range, &fftw_plans, &fftw_arr, &aux_arr);

     return(EXIT_SUCCESS);

}


// signal handler to save state and exit before end
static void sig_handler(int signo)
{
     if (signo == SIGTERM || signo == SIGUSR1) save_state = 1;
}
