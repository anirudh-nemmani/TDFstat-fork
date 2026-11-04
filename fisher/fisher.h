#ifndef __FISHER_H__
#define __FISHER_H__

#include "struct.h"

// Parameters, in ParamIndex order: FREQ, FDOT, RA, DEC, IOTA, PSI, PHASE, H0
#define NPARAMS 8

// Fisher matrix of the detector network
void fisher_matrix( Fisher_settings *sett,
                    Command_line_opts *opts,
                    Aux_arrays *aux_arr,
                    Signal_params *sgnl_params);

// Systematic biases for each detector, from dA and dPsi in opts->sysfile
void systematics( Fisher_settings *sett,
                  Command_line_opts *opts,
                  Aux_arrays *aux_arr,
                  Signal_params *sgnl_params);

#endif
