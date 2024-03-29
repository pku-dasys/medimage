# Makefile template by Guojie Luo

ROOT_BIN = ./bin
ROOT_SRC = ./src
ROOT_LIB = ./lib

ROOT_DIR = .

ROOT_VPATH = $(foreach DIR, $(ROOT_DIR), $(ROOT_SRC)/$(DIR))
ROOT_INCL = $(foreach DIR, $(ROOT_DIR), -I$(ROOT_SRC)/$(DIR))

ROOT_CCFLAGS = $(ROOT_INCL) -I$(ROOT_LIB) -std=c89 -D_GNU_SOURCE\
							-Wall -m64 -march=native -funroll-all-loops -fprefetch-loop-arrays -ffast-math -fopenmp

ROOT_LDFLAGS = -lm -lrt \
							-m64 -march=native -funroll-all-loops -fprefetch-loop-arrays -ffast-math -fopenmp

ifeq ($(build), debug)
ROOT_POSTFIX=-debug
ROOT_OPT = -g
else
ifeq ($(build), gprof)
ROOT_POSTFIX=-gprof
ROOT_OPT = -pg
else
ROOT_POSTFIX=
ROOT_OPT = -O3
endif
endif

ROOT_CC = gcc ${ROOT_OPT}
ROOT_LNK = gcc ${ROOT_OPT}

ROOT_NAME = MSbeam
ROOT_BUILD = ./build$(ROOT_POSTFIX)

all ::
	if [ ! -d $(ROOT_BUILD) ]; then mkdir $(ROOT_BUILD) || exit 1 ; fi
	if [ ! -d $(ROOT_BIN) ]; then mkdir $(ROOT_BIN) || exit 1 ; fi

vpath %.c $(ROOT_VPATH)

include $(foreach TOOL, $(ROOT_NAME), $(ROOT_SRC)/$(TOOL).mk)

$(ROOT_BUILD)/%.o :: %.c
	if $(ROOT_CC) -c $(ROOT_CCFLAGS) -MD -MT $@ -MP -MF $(ROOT_BUILD)/$*.CXXd $< -o $@; \
		then mv -f "$(ROOT_BUILD)/$*.CXXd" "$(ROOT_BUILD)/$*.d" ; \
		else rm -f "$(ROOT_BUILD)/$*.CXXd" ; exit 1 ; fi
