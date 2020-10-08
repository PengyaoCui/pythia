SHELL = /bin/bash

DMODE = 0 # 0: TEST
          # 1: CERN
          # 2: CCNU

SRC := $(wildcard *.C)
TAG := $(patsubst %.C, %, $(SRC))
INC := $(wildcard inc/*.h)

ROOT_CXXFLAGS := `root-config --cflags`
ROOT_LDLIBS   := `root-config --libs` -lEG
ROOT_LDFLAGS  := `root-config --ldflags`

PYTHIA8_CXXFLAGS := `pythia8-config --cxxflags`
PYTHIA8_LIBDIR   := `pythia8-config --libdir`
PYTHIA8_LDLIBS   := -L$(PYTHIA8_LIBDIR) -Wl,-rpath $(PYTHIA8_LIBDIR) -lpythia8

FASTJET_CXXFLAGS := `fastjet-config --cxxflags`
FASTJET_LDLIBS   := `fastjet-config --libs`

EXTRA_CXXFLAGS := -Iinc -D CMODE=$(DMODE)

CXX      := `root-config --ld`
CXXFLAGS := $(ROOT_CXXFLAGS) $(PYTHIA8_CXXFLAGS) $(FASTJET_CXXFLAGS) $(EXTRA_CXXFLAGS)
LDLIBS   := $(ROOT_LDLIBS)   $(PYTHIA8_LDLIBS)   $(FASTJET_LDLIBS)

LDFLAGS := -g -O3 -W -Wall -Wshadow -fPIC $(ROOT_LDFLAGS)

all : $(TAG)
.PHONY : all

% : %.C $(INC)

clean :
	@$(RM) -r $(TAG) *.dSYM
.PHONY : clean

test-clean :
	@$(RM) *.root
.PHONY : test-clean

dist-clean : clean test-clean
.PHONY : dist-clean
