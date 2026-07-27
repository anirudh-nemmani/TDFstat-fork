# arch to compile for; min. required for SLEEF is haswell (AVX2+FMA)
ARCH = cascadelake

# SIMD = [AVX2, AVX512, NONE] ; explicit optimizations with SIMD instructions.
# SIMD=NONE implies usege of sincos function, otherwise the SLEEF library is used.
# It does not affect SIMD optimizations by the compiler and FFTW3 dispatchers.
SIMD = AVX2
# legacy option: USE_AVX2_SINCOS

# PREC = [FLOAT, DOUBLE] precission used in the spindown loop and for output
PREC = FLOAT

#FFTW_DIR = $(CONDA_PREFIX)
FFTW_DIR = /opt/fftw/3.3.10-gcc11-mvapich2/

#-------------------- EDIT BELOW WITH CAUTION --------------

INIPARSER_DIR = ../../../utils/iniparser

ifeq ($(strip $(SIMD)), NONE)
	PERF = sincos-${PREC}
	SINCOS = GNUSINCOS
	VCFLAGS =
	VLDFLAGS =
	VLIBS =
else
	PERF = ${SIMD}-${PREC}
	SINCOS = SLEEF
	SPATH = lib/sleef-3.6
	VCFLAGS  = -DUSE_${SIMD} -I$(SPATH)/include
	VLDFLAGS = -L$(SPATH)/lib64
	VLIBS = -lsleef
endif

ifeq ($(strip $(PREC)), FLOAT)
	FFTW_FLOAT_LIBS = -lfftw3f -lfftw3f_omp
else
	FFTW_FLOAT_LIBS =
endif


#TARGET = gwsearch-${ARCH}-$(shell echo $(PERF) | tr A-Z a-z)-oldtrig
TARGET = gwsearch-${ARCH}-$(shell echo $(PERF) | tr A-Z a-z)

VER := "\"$(shell git rev-parse HEAD)\""

CC = gcc -std=gnu17 -static-libgcc -fopenmp #-flto
#-flto ?? speeds up a little for AVX512 but slows down for AVX2

# TIMERS=1 measure total time, TIMERS=3 measure spindown loop time
# PREFIX is the default output directory
CFLAGS = -DPREFIX="./candidates" -DCODEVER=$(VER) -DTIMERS=3 -D$(SINCOS) -UUSE_LOCKING  -DCOMP_${PREC} -DSCI_RUN=${SCI_RUN} ${VCFLAGS} #-DVERBOSE
CFLAGS += -I$(FFTW_DIR)/include -I$(INIPARSER_DIR)/src
CFLAGS += -g -Wall -Wno-unused \
	  -O3 -ffast-math  -funsafe-loop-optimizations -funroll-loops \
	  -march=$(ARCH) -mtune=$(ARCH)


LDFLAGS = -L$(FFTW_DIR)/lib ${VLDFLAGS} -L$(INIPARSER_DIR)
LDLIBS = -Wl,--dynamic-linker=/lib64/ld-linux-x86-64.so.2 \
	 -static -lfftw3 -lfftw3_omp $(FFTW_FLOAT_LIBS) ${VLIBS} \
	 -lgsl -lgslcblas -liniparser \
	 -Wl,-Bdynamic -lc -lm -lrt -lhdf5 -lhdf5_hl



all: veclib $(TARGET)

ifeq ($(strip $(SINCOS)),SLEEF)
veclib:
	@echo -e "\nUsing prebuild SLEEF library from: $(SPATH)\n"
endif

SRC = main.c jobcore.c auxi.c settings.c init.c  timer.c spinmod.c io.c common.c
INC = jobcore.h settings.h auxi.h timer.h init.h struct.h io.h common.h
DEPS = $(SRC) $(INC) Makefile

$(TARGET): $(DEPS)
	$(CC) $(LDFLAGS) -o $@ $(SRC) $(CFLAGS) $(LDLIBS)

clean:
	rm -f *.o

uninstall: clean
	rm -f $(TARGET)

.PHONY: all $(TARGET) veclib clean uninstall

#----------------------------------------------

sigen: sigen.c auxi.c settings.c
	$(CC) -o $@ $^ -Wl,-Bdynamic -lc -lrt -lm -lgsl -lgslcblas

#----------------------------------------------

# Standalone tool: reads a search .ini file and writes out vlines.dat
# (veto lines found in the given band). Links the same init.c and settings.c
# as the search itself, so no search code is duplicated.
#
# vlinegen.c does keep its own read_lines(), line_in_band(), narrow_down_band()
# and lines_veto_fraction(). Those four names also exist in settings.c, so that
# file is compiled here with the names renamed, which leaves the copies in
# vlinegen.c as the only definitions and avoids a clash at link time.
VLG_RENAME = -Dread_lines=settings_read_lines \
	     -Dline_in_band=settings_line_in_band \
	     -Dnarrow_down_band=settings_narrow_down_band \
	     -Dlines_veto_fraction=settings_lines_veto_fraction

VLG_CFLAGS = -DCOMP_FLOAT -I$(FFTW_DIR)/include -I$(INIPARSER_DIR)/src
VLG_LIBS = -L$(FFTW_DIR)/lib -L$(INIPARSER_DIR) \
	   -Wl,-Bdynamic -lc -liniparser -lgsl -lgslcblas \
	   -lfftw3 -lfftw3_omp -lfftw3f -lfftw3f_omp -lhdf5 -lhdf5_hl -lm -lrt

vlinegen: vlinegen.c auxi.c init.c settings.c common.c \
	  settings.h struct.h auxi.h init.h common.h Makefile
	$(CC) $(VLG_CFLAGS) $(VLG_RENAME) -c -o settings-vlg.o settings.c
	$(CC) $(VLG_CFLAGS) -o $@ vlinegen.c auxi.c init.c common.c settings-vlg.o \
	    $(VLG_LIBS)
	rm -f settings-vlg.o
