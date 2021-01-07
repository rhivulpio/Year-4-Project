Target  = Process
ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --libs)

DELPHESEXT = $(DELPHES)/external
DELPHESLIB = $(DELPHES)/lib


all:$(Target)

Process: Process.cxx Process.h  $(Objects)
	g++ -g -o $@ Process.cxx $(Objects) -L$(DELPHESLIB) -I$(DELPHES) -I$(DELPHESEXT) $(ROOTFLAGS) $(ROOTLIBS) $(DELPHES)/libDelphes.so

clean:
	rm -f Process *.o
