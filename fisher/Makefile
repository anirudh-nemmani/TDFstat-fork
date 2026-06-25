# arch to compile for; min. required for SLEEF is haswell (AVX2+FMA)
# ARCH = sapphirerapids
ARCH = cascadelake

# SIMD = [AVX2, AVX512, NONE] ; explicit optimizations with SIMD instructions.
# SIMD=NONE implies usege of sincos function, otherwise the SLEEF library is used.
# It does not affect SIMD optimizations by the compiler and FFTW3 dispatchers.
SIMD = AVX2
# legacy option: USE_AVX2_SINCOS

# PREC = [FLOAT, DOUBLE] precission used in the spindown loop and for output
PREC = FLOAT

#-------------------- EDIT BELOW WITH CAUTION --------------

INIPARSER_DIR = ../utils/iniparser

ifeq ($(strip $(SIMD)), NONE)
	PERF = sincos-${PREC}
	SINCOS = GNUSINCOS
	VCFLAGS =
	VLDFLAGS =
	VLIBS =
else
	PERF = ${SIMD}-${PREC}
	SINCOS = SLEEF
	SPATH = ../search/network/openmp/lib/sleef-3.6
	VCFLAGS  = -DUSE_${SIMD} -I$(SPATH)/include
	VLDFLAGS = -L$(SPATH)/lib64
	VLIBS = -lsleef
endif


#TARGET = gwsearch-${ARCH}-$(shell echo $(PERF) | tr A-Z a-z)-oldtrig
TARGET = gwfisher-${ARCH}-$(shell echo $(PERF) | tr A-Z a-z)

VER := "\"$(shell git rev-parse HEAD)\""

CC = gcc -std=gnu17 -static-libgcc -fopenmp #-flto
#-flto ?? speeds up a little for AVX512 but slows down for AVX2

# TIMERS=1 measure total time, TIMERS=3 measure spindown loop time
# PREFIX is the default output directory
CFLAGS = -DPREFIX="./candidates" -DCODEVER=$(VER) -DTIMERS=3 -D$(SINCOS) -UUSE_LOCKING  -DCOMP_${PREC} -DSCI_RUN=${SCI_RUN} ${VCFLAGS} #-DVERBOSE
CFLAGS += -I$(INIPARSER_DIR)/src
CFLAGS += -g -Wall -Wno-unused \
	  -O3 -ffast-math  -funsafe-loop-optimizations -funroll-loops \
	  -march=$(ARCH) -mtune=$(ARCH)


LDFLAGS = ${VLDFLAGS} -L$(INIPARSER_DIR)
LDLIBS = -Wl,--dynamic-linker=/lib64/ld-linux-x86-64.so.2 \
	 -static \
	 -lgsl -lgslcblas -liniparser \
	 -Wl,-Bdynamic -lc -lm -lrt -lhdf5 -lhdf5_hl



all: veclib $(TARGET)

ifeq ($(strip $(SINCOS)),SLEEF)
veclib:
	@echo -e "\nUsing prebuild SLEEF library from: $(SPATH)\n"
endif

SRC = main.c fisher.c auxi.c settings.c init.c
INC = fisher.h settings.h auxi.h init.h struct.h
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
