# Need this to get SHAREDSUFFIX (e.g. dylib or so)
#-include ../config.mk

EX = cHNLdecay
CC = gcc
CXX = g++
RM = rm -f
CPPFLAGS = -g $(shell root-config --cflags) -I/usr/include
LDFLAGS = -g $(shell root-config --ldflags)
LDLIBS = $(shell root-config --libs) -lmpfr -lgmp -L/usr/lib -L/usr/lib64
LDLIBS_DBG = $(shell root-config --libs) -lmpfr -lgmp -L/usr/lib -L/usr/lib64 -lprofiler -ltcmalloc

SRCS = $(wildcard *.cxx)
#SRCS = auxfunctions.cxx HNL.cxx Logger.cxx partialWidths.cxx plots.cxx prodFromBmesons.cxx
OBJS = $(SRCS:.cxx=.o)


#all: $(EX)

all: $(EX)

debug: $(OBJS)
	$(CXX) $(LDFLAGS) -o cHNLdecay $(EX) $(OBJS) $(LDLIBS_DBG) -ggdb

cHNLdecay: $(OBJS) cHNLdecay.o
	$(CXX) $(LDFLAGS) -o cHNLdecay $(OBJS) $(LDLIBS)
#cHNLproduction: $(OBJS) cHNLproduction.o
#	$(CXX) $(LDFLAGS) -o cHNLproduction $(OBJS) $(LDLIBS)
	
	

# To obtain object files
%.o: %.cxx
	$(CXX) -c $(CPPFLAGS) $< -o $@

# To remove generated files
clean:
	rm -f $(EXEC) $(OBJS) cHNLdecay.o cHNLproduction.o
