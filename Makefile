CC=g++
OPTFLAGS=-O2 -march=native
CFLAGS=-Wall --std=c++11 -fpermissive $(OPTFLAGS)
INCLUDES = -I.
LFLAGS= -lssl -lcrypto

# OPTIONS
# change this for different optimizations
# This tells which folder to take scalarmult and radix from
OPTDIR =SSE
# OPTDIR =StdC-optimize
# OPTDIR =straight
# This tells benchmark which include to do

ifeq ($(OPTDIR),SSE)
OPTOPT =SSE
endif
ifeq ($(OPTDIR),StdC-optimize)
OPTOPT =STDC
endif
ifeq ($(OPTDIR),straight)
OPTOPT =BASE
endif

RADIXOPT =R25

BUILDSUBDIR =build

# directories
BUILDDIR =$(OPTDIR)/$(BUILDSUBDIR)
RADIXDIR =$(OPTDIR)/radix
RADIX25DIR =$(OPTDIR)/radix25
SCALARMULDIR =$(OPTDIR)/scalarmult

# SOURCES
RADIXSRCS = $(addprefix $(RADIXDIR)/, radix51_common.cpp radix51_addsub.cpp radix51_mul.cpp radix51_shift.cpp radix51_sqr.cpp radix51_inv.cpp radix51.cpp radix51_expmod.cpp )
RADIX25SRCS = $(addprefix $(RADIX25DIR)/, radix25_common.cpp radix25_addsub.cpp radix25_mul.cpp radix25_shift.cpp radix25_sqr.cpp radix25_inv.cpp radix25.cpp radix25_expmod.cpp )
SCALARMULTSRCS = $(addprefix $(SCALARMULDIR)/, scalarmult.cpp)
SRCS = $(RADIXSRCS) $(RADIX25SRCS) $(SCALARMULTSRCS) converter.cpp util.cpp ed25519.cpp main.cpp
BENCHSRCS = $(RADIXSRCS) $(RADIX25SRCS) $(SCALARMULTSRCS) converter.cpp util.cpp ed25519.cpp main.cpp benchmark.cpp
# OBJECTS
OBJS = radix51.o scalarmult.o ed25519.o main.o util.o
OBJS25 = radix25.o scalarmult25.o ed25519_25.o main.o util.o
BENCHOBJS = radix51.o scalarmult.o util.o ed25519.o benchmark.o
BENCH25OBJS = radix25.o scalarmult25.o util.o ed25519_25.o benchmark25.o
# BINARY
MAIN = proj
MAIN25 = proj25

.PHONY: depend clean

all:    $(MAIN) $(MAIN25) benchmark bench25

info:
	@echo ""
	@echo "** Options **"
	@echo "OPTDIR=$(value OPTDIR)	|	OPTOPT=$(value OPTOPT)	|	RADIXOPT=$(value RADIXOPT)"
	@echo "- all:	proj proj25 benchmark bench25"
	@echo "- proj:		builds main file with current options and radix51."
	@echo "- proj25:	builds main file with current options and radix25."
	@echo "- benchmark:	builds benchmark with current options and radix51."
	@echo "- bench25:	builds benchmark with current options and using radix25."
	@echo "*********"

$(MAIN): $(OBJS)
	$(CC) $(CFLAGS) -D $(OPTOPT) $(INCLUDES) -o $(BINDIR)$(MAIN) $(addprefix $(BUILDDIR)/, $(OBJS)) $(LFLAGS)

$(MAIN25): $(OBJS25)
	$(CC) $(CFLAGS) -D $(OPTOPT) -D $(RADIXOPT) $(INCLUDES) -o $(BINDIR)$(MAIN25) $(addprefix $(BUILDDIR)/, $(OBJS25)) $(LFLAGS)

benchmark:$(BENCHOBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -D $(OPTOPT) -o $(BINDIR)benchmark $(addprefix $(BUILDDIR)/, $(BENCHOBJS)) $(LFLAGS)

bench25:$(BENCH25OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -D $(OPTOPT) -D $(RADIXOPT) -o $(BINDIR)bench25 $(addprefix $(BUILDDIR)/, $(BENCH25OBJS)) $(LFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) -D $(OPTOPT) $(INCLUDES) -c $<  -o $(addprefix $(BUILDDIR)/, $@) $(LFLAGS)

benchmark25.o: benchmark.cpp
	$(CC) $(CFLAGS) -D $(OPTOPT) -D $(RADIXOPT) $(INCLUDES) -c $<  -o $(addprefix $(BUILDDIR)/, $@) $(LFLAGS)

ed25519.o: ed25519.cpp
	$(CC) $(CFLAGS) -D $(OPTOPT) $(INCLUDES) -c $<  -o $(addprefix $(BUILDDIR)/, $@) $(LFLAGS)

ed25519_25.o: ed25519.cpp
	$(CC) $(CFLAGS) -D $(OPTOPT) -D $(RADIXOPT) $(INCLUDES) -c $<  -o $(addprefix $(BUILDDIR)/, $@) $(LFLAGS)

radix51.o:
	$(MAKE) radix51.o -C $(OPTDIR) OPTFLAGS="$(OPTFLAGS)" BUILDSUBDIR=$(BUILDSUBDIR)

radix25.o:
	$(MAKE) radix25.o -C $(OPTDIR) OPTFLAGS="$(OPTFLAGS)" BUILDSUBDIR=$(BUILDSUBDIR)

scalarmult.o:
	$(MAKE) scalarmult.o -C $(OPTDIR) OPTFLAGS="$(OPTFLAGS)" BUILDSUBDIR=$(BUILDSUBDIR)

scalarmult25.o:
	$(MAKE) scalarmult25.o -C $(OPTDIR) OPTFLAGS="$(OPTFLAGS)" BUILDSUBDIR=$(BUILDSUBDIR)

clean:
	$(RM) $(BUILDDIR)/*.o *~ $(BINDIR)$(MAIN) $(BINDIR)$(MAIN25) $(BINDIR)benchmark $(BINDIR)bench25

depend: $(BENCHSRCS)
	makedepend $(INCLUDES) $^
# DO NOT DELETE THIS LINE -- make depend needs it
