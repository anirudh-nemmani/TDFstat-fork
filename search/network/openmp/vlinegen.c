#include <math.h>       // M_PI
#include <stdio.h>      // printf, perror
#include <unistd.h>

#include <errno.h>      // errno, ENOENT
#include <sys/stat.h>   // stat, mkdir

#include "common.h"     // read_detssb, narrow_down_band, extract_band_vlines
#include "struct.h"     // Search_settings, Command_line_opts, Detector_settings
#include "settings.h"   // detectors_settings, search_settings
#include "init.h"       // read_ini_file, read_grid

/*
 * vlinegen - reads a search .ini file and writes vlines_${BAND}.dat, the list of
 * instrumental veto lines found in the given band, plus the veto fraction.
 *
 * main() follows the same initialization sequence as the search (main.c).
 */

/* Global array of detectors (network); declared extern in struct.h,
 * defined once here since vlinegen does not link against main.c */
Detector_settings ifo[MAX_DETECTORS];

int main(int argc, char *argv[])
{

    Search_settings sett;
    Command_line_opts opts;
    struct stat buffer;
    int i;

    // Command line options (reads the .ini file given as argv[1])
    read_ini_file(&sett, &opts, argc, argv);

    // prevent extract_band_vlines from applying lines, we only generate them.
    opts.veto_flag = 0;

    // Output directory handling (create it if missing, like main.c does)
    if (stat(opts.outdir, &buffer) == -1) {
        if (errno == ENOENT) {
            if (mkdir(opts.outdir, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH) == -1) {
                perror(opts.outdir);
                return 1;
            }
        } else {
            perror(opts.outdir);
            return 1;
        }
    }

    // Detector network settings
    detectors_settings(&sett, &opts);

    // Grid data (sets sett.fftpad, needed by search_settings())
    read_grid(&sett, &opts);

    // Search settings (sets sett.N, sett.B, sett.numlines_band, ...
    // all of them used by read_lines())
    search_settings(&sett);

    // Reading the ephemerids, i.e. ifo[].sig.DetSSB needed by read_lines().
    // Only this part of init_arrays() is used here; the input time series and
    // the search work arrays are not needed to generate the veto lines.
    for (i=0; i<sett.nifo; i++)
        read_detssb(&sett, &opts, i);

    // Narrowing-down the band (excludes the edges according to
    // opts.narrowdown parameter); adds two lines
    if (opts.narrowdown < 0.5*M_PI) narrow_down_band(&sett, &opts);

    char band_vl_file[FNAME_LENGTH];
    sprintf(band_vl_file, "%s/vlines_%04d%s.dat", opts.outdir, opts.band, opts.label);
    if (access(band_vl_file, R_OK) == 0) {
        printf("[Warning] Overwriting file %s !\n", band_vl_file);
    }
    // Reading veto lines data from external files, writes vlines##.dat
    printf("Reading veto files...\n");
    extract_band_vlines(&sett, &opts, band_vl_file);

    return 0;
}
