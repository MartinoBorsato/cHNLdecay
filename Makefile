# Need this to get SHAREDSUFFIX (e.g. dylib or so)
#-include ../config.mk

EX = cHNLdecay
CPPFLAGS += -g $(shell root-config --cflags)
LDFLAGS += -g $(shell root-config --ldflags)
LDLIBS += $(shell root-config --libs) -lmpfr -lgmp
LDLIBS_DBG += $(shell root-config --libs) -lmpfr -lgmp -lprofiler -ltcmalloc

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
