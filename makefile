CXX = $(shell root-config --cxx)
CPPFLAGS = -isystem$(shell root-config --incdir) -I inc
CXXFLAGS = -Wall -Wextra -pedantic -O2 -Wshadow -fPIC $(shell root-config --cflags)
LD = $(shell root-config --ld)
LDFLAGS = $(shell root-config --ldflags)
LDLIBS =  $(shell root-config --glibs) -lMinuit2

VPATH = inc:src:obj

OBJECTS = obj/METSigFit.o

COMPILE = $(CXX) $(CXXFLAGS) $(CPPFLAGS) -c
LINK = $(LD) $(LDFLAGS)
LINKEND = $(OBJECTS) $(LDLIBS)

# add your executable name here
all : DoFit

DoFit : DoFit.o $(OBJECTS)
	$(LINK) -o DoFit DoFit.o $(LINKEND)
DoFit.o: DoFit.C METSigFit.h
	$(COMPILE) DoFit.C

clean:
	-rm -f DoFit obj/*.o *.o

obj/METSigFit.o : METSigFit.C METSigFit.h
	$(COMPILE) src/METSigFit.C -o obj/METSigFit.o

.PHONY : clean
