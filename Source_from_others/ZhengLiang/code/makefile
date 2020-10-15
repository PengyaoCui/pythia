OBJS = AnaPythia.o 
EXE = AnaPythia.exe

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)

#PYTHIAINC=${HOME}/lib/pythia8/include/
#PYTHIALIBS=${HOME}/lib/pythia8/lib/
PYTHIALIBS=${HOME}/alice/sw/ubuntu1804_x86-64/pythia/alice-1/lib/
PYTHIAINC=${HOME}/alice/sw/ubuntu1804_x86-64/pythia/alice-1/include/

INCFLAGS = -I$(ROOTSYS)/include
LDFLAGS = -L$(ROOTSYS)/lib

CXX = g++

all: $(EXE)

$(EXE): $(OBJS)
	$(CXX) -I$(PYTHIAINC) -I$(FASTJETINC) $(ROOTCFLAGS) -o $(EXE) $(OBJS) $(ROOTLIBS) -L$(PYTHIALIBS)/ -lpythia8 

AnaPythia.o: AnaPythia.C
	$(CXX) -g -I$(PYTHIAINC)  $(ROOTCFLAGS) -c AnaPythia.C


.PHONY : clean
clean :
	-rm -f $(OBJS) $(EXE)
