/*
Compile using:

gcc -std=gnu17 -static-libgcc -fopenmp  -L/opt/fftw/3.3.10-gcc11-mvapich2//lib -L/opt/sleef/3.6/lib64 -L../../../utils/iniparser -o write-signal write-signal.c jobcore.c auxi.c settings.c init.c  timer.c spinmod.c io.c -DPREFIX="./candidates" -DCODEVER="\"76de1ff38fb4c83399654aadecd944a1bad80309\"" -DTIMERS=1 -DSLEEF -UUSE_LOCKING  -DCOMP_FLOAT -DSCI_RUN= -DUSE_AVX2 -I/opt/sleef/3.6/include  -I/opt/fftw/3.3.10-gcc11-mvapich2//include -I../../../utils/iniparser/src -g -Wall -Wno-unused -O3 -ffast-math  -funsafe-loop-optimizations -funroll-loops -march=sapphirerapids -mtune=sapphirerapids -Wl,--dynamic-linker=/lib64/ld-linux-x86-64.so.2 -static -lfftw3 -lfftw3_omp -lfftw3f -lfftw3f_omp -lsleef -lgsl -lgslcblas -liniparser -Wl,-Bdynamic -lc -lm -lrt -lhdf5 -lhdf5_hl

*/

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

     // Reading veto lines data from external files
     printf("Reading veto files...\n");
     read_lines(&sett, &opts);
     if (opts.gen_vlines_flag) exit(EXIT_SUCCESS);

     // Amplitude modulation functions for each detector
     for(i=0; i<sett.nifo; i++)
          rogcvir(&ifo[i]);

    // writing the signal into the memory
    add_signal(&sett, &opts, &aux_arr);

    // Writing ifo[n].sig.xDat into a file
    for (i = 0; i < sett.nifo; i++) {
        FILE *file = fopen(ifo[i].xdatname, "wb");
        if (file == NULL) {
                perror("Error opening file for writing");
                exit(EXIT_FAILURE);
        }

        size_t data_size = sett.N * sizeof(float);
        if (fwrite(ifo[i].sig.xDat, 1, data_size, file) != data_size) {
                perror("Error writing data to file");
                fclose(file);
                exit(EXIT_FAILURE);
        }
        fclose(file);
    }

}

// signal handler to save state and exit before end
static void sig_handler(int signo)
{
     if (signo == SIGTERM || signo == SIGUSR1) save_state = 1;
}
