CC=g++
OPTFLAGS=-O3 -march=native
CFLAGS=-Wall --std=c++11 -fpermissive $(OPTFLAGS)
INCLUDES = -I.
LFLAGS= -lssl -lcrypto

BUILDSUBDIR = build

BUILDDIR = $(BUILDSUBDIR)
RADIXDIR = radix
RADIX25DIR = radix25

# SOURCES
RADIXSRCS = $(addprefix $(RADIXDIR)/, radix51_common.cpp radix51_addsub.cpp radix51_mul.cpp radix51_shift.cpp radix51_sqr.cpp radix51_inv.cpp radix51.cpp radix51_expmod.cpp )
RADIX25SRCS = $(addprefix $(RADIX25DIR)/, radix25_common.cpp radix25_addsub.cpp radix25_mul.cpp radix25_shift.cpp radix25_sqr.cpp radix25_inv.cpp radix25.cpp radix25_expmod.cpp )
SCALARMULTSRCS = scalarmult/scalarmult.cpp
TESTSRC = testsrc/test.cpp
TESTS25SRC = testsrc/testsradix25.cpp
# OBJECTS
TESTOBJS = radix51.o scalarmult.o tests.o
TEST25OBJS = radix25.o scalarmult25.o testsradix25.o
CONVERTEROBJS = radix51.o converter.o
# BINARY

.PHONY: depend clean

all:    test test25

test:	$(TESTOBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o test $(addprefix $(BUILDDIR)/, $(TESTOBJS)) $(LFLAGS)

test25:	$(TEST25OBJS)
	$(CC) $(CFLAGS) -D R25 $(INCLUDES) -o test25 $(addprefix $(BUILDDIR)/, $(TEST25OBJS)) $(LFLAGS)

converter:	$(CONVERTEROBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o converter $(addprefix $(BUILDDIR)/, $(CONVERTEROBJS))

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $(addprefix $(BUILDDIR)/, $@) $(LFLAGS)

radix51.o: $(RADIXSRCS)
	$(CC) $(CFLAGS) $(INCLUDES) -c radix/unity.cpp -o $(BUILDDIR)/$@ $(LFLAGS)

radix25.o: $(RADIX25SRCS)
	$(CC) $(CFLAGS) -D R25 $(INCLUDES) -c radix25/unity.cpp -o $(BUILDDIR)/$@ $(LFLAGS)

scalarmult.o: $(SCALARMULTSRCS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $(BUILDDIR)/$@ $(LFLAGS)

scalarmult25.o: $(SCALARMULTSRCS)
	$(CC) $(CFLAGS) -D R25 $(INCLUDES) -c $<  -o $(BUILDDIR)/$@ $(LFLAGS)

tests.o: testsrc/tests.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $(BUILDDIR)/$@ $(LFLAGS)

testsradix25.o: testsrc/testsradix25.cpp
	$(CC) $(CFLAGS) -D R25 $(INCLUDES) -c $< -o $(BUILDDIR)/$@ $(LFLAGS)

converter.o: testsrc/converter.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $(BUILDDIR)/$@ $(LFLAGS)

clean:
	$(RM) $(BUILDDIR)/*.o *~ test converter test25

depend: $(BENCHSRCS)
	makedepend $(INCLUDES) $^
# DO NOT DELETE THIS LINE -- make depend needs it
