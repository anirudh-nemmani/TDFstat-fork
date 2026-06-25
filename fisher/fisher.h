#ifndef __FISHER_H__
#define __FISHER_H__

#include "struct.h"

void fisher_matrix(
                    Fisher_settings *sett,
                    Command_line_opts *opts,
                    Aux_arrays *aux_arr,
                    Signal_params *sgnl_params);

int detector_settings_deep_copy(Detector_settings *dst,
                                const Detector_settings *src,
                                Fisher_settings *sett);

int detector_array_deep_copy(Detector_settings *dst_arr,
                             const Detector_settings *src_arr,
                             Fisher_settings *sett);

void detector_array_free_owned(Detector_settings *arr,
                               Fisher_settings *sett);

void perturb_params(Signal_params *out,
                    Signal_params *base,
                    ParamIndex p,
                    double dp);

#endif
