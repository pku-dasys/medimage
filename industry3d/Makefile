# Makefile template by Guojie Luo

ROOT_BIN = ./bin
ROOT_SRC = ./src
ROOT_LIB = ./lib

ROOT_DIR = . testcase

ROOT_VPATH = $(foreach DIR, $(ROOT_DIR), $(ROOT_SRC)/$(DIR))
ROOT_INCL = $(foreach DIR, $(ROOT_DIR), -I$(ROOT_SRC)/$(DIR))

ROOT_CCFLAGS = $(ROOT_INCL) -I$(ROOT_LIB) \
							 -Wreorder -Wformat -Wunused \
							 -I./include -fopenmp -std=c++11

ROOT_LDFLAGS = -L$(ROOT_LIB) -fopenmp -lboost_system -lboost_filesystem

ifeq ($(build), debug)
ROOT_POSTFIX=.debug
ROOT_OPT = -g
else
ifeq ($(build), prof)
ROOT_POSTFIX=.prof
ROOT_OPT = -pg -g
else
ROOT_POSTFIX=.release
ROOT_OPT = -O3
endif
endif

ROOT_NAME = ct3d generate
ROOT_BUILD = ./build$(ROOT_POSTFIX)

all ::
	if [ ! -d $(ROOT_BUILD) ]; then mkdir $(ROOT_BUILD) || exit 1 ; fi
	if [ ! -d $(ROOT_BIN) ]; then mkdir $(ROOT_BIN) || exit 1 ; fi
	mkdir -p $(foreach DIR, $(ROOT_DIR), $(ROOT_BUILD)/$(DIR))

ROOT_CC = g++ ${ROOT_OPT}
ROOT_LNK = g++ ${ROOT_OPT}

vpath %.cpp $(ROOT_VPATH)

include $(foreach TOOL, $(ROOT_NAME), $(ROOT_SRC)/$(TOOL).mk)

$(ROOT_BUILD)/%.o :: %.cpp
	if $(ROOT_CC) -c $(ROOT_CCFLAGS) -MD -MT $@ -MP -MF $(ROOT_BUILD)/$*.CXXd $< -o $@; \
		then mv -f "$(ROOT_BUILD)/$*.CXXd" "$(ROOT_BUILD)/$*.d" ; \
		else rm -f "$(ROOT_BUILD)/$*.CXXd" ; exit 1 ; fi
