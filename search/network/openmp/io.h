//#include <stdlib.h>
//#include <math.h>
//#include <complex.h>

#include <complex.h>
#include "auxi.h"
#include "struct.h"

int trig_h5_init (char *outname, Command_line_opts *opts, Search_settings *sett,
     Search_range *s_range, Trigger *sgnlv);

int trig_h5_extend(char *outname, Trigger *sgnlv, int sgnlv_size);
