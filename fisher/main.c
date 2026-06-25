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
//#include <getopt.h>
#include <gsl/gsl_linalg.h>
#include <time.h>
#include <dirent.h>
#include <signal.h>

#include "auxi.h"
#include "struct.h"
#include "settings.h"
#include "fisher.h"
#include "init.h"

#ifndef CODEVER
#define CODEVER unknown
#endif

Detector_settings ifo[MAX_DETECTORS];
volatile sig_atomic_t save_state = 0;

int main (int argc, char* argv[])
{

     Command_line_opts opts;
     Fisher_settings sett;
     Aux_arrays aux_arr;
     Signal_params sgnl_params;
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

     // Fisher settings
     fisher_settings(&sett);

     // Array initialization and reading the ephemerids
     init_arrays(&sett, &opts, &aux_arr);

     // Amplitude modulation functions for each detector
     for (i=0; i<sett.nifo; i++)
          rogcvir(&ifo[i]);

     // Generate signal and store it
     read_signal_file(&sgnl_params, &opts);
     // signal_gen(&sett, &opts, &aux_arr, &sgnl_params, ifo);

     // Fisher matrix calculation
     fisher_matrix(&sett, &opts, &aux_arr, &sgnl_params);
}

// signal handler to save state and exit before end
static void sig_handler(int signo)
{
     if (signo == SIGTERM || signo == SIGUSR1) save_state = 1;
}
