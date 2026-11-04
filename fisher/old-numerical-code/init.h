#ifndef __INIT_H__
#define __INIT_H__

#include "struct.h"

void read_ini_file(
                    Fisher_settings *sett,
                    Command_line_opts *opts,
                    int argc,
                    char* argv[]);

void init_arrays(
                    Fisher_settings *sett,
                    Command_line_opts *opts,
                    Aux_arrays *aux_arr);

void read_signal_file(
                    Signal_params *sgnl_params,
                    Command_line_opts *opts);

void signal_gen(
                    Fisher_settings *sett,
                    Command_line_opts *opts,
                    Aux_arrays *aux_arr,
                    Signal_params *sgnl_params,
                    Detector_settings *ifo);

void plan_fftw(
                    Fisher_settings *sett,
                    Command_line_opts *opts,
                    Aux_arrays *aux_arr);

void cleanup(
                    Fisher_settings *sett,
                    Command_line_opts *opts,
                    Aux_arrays *aux);

#endif
