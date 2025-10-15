//#include <stdlib.h>
//#include <math.h>
//#include <complex.h>

#include "auxi.h"
#include "struct.h"

int hdfout_init (char *outname, Command_line_opts *opts, Search_settings *sett,
     Search_range *s_range, Trigger *sgnlv);

int hdfout_extend(char *outname, Trigger *sgnlv, int sgnlv_size);
int hdfout_finalize(char *outname, int totsgnl, double time_elapsed, int nthreads);
