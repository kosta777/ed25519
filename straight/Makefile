CC=g++
CFLAGS=-Wall --std=c++11 -O3 -fno-tree-vectorize
INCLUDES = -I.
LFLAGS= -lssl -lcrypto

BUILDDIR = build

SRCS = radix51.cpp scalarmult.cpp converter.cpp util.cpp ed25519.cpp main.cpp
TESTSRC = testsrc/test.cpp
RADIXSRC = radix51.cpp
# BENCHSRCS = radix51.cpp scalarmult.cpp converter.cpp util.cpp ed25519.cpp main.cpp benchmark.cpp
# OBJS = radix51.o scalarmult.o ed25519.o main.o util.o
TESTOBJS = radix51.o scalarmult.o tests.o
# BENCHOBJS = radix51.o scalarmult.o util.o ed25519.o benchmark.o

CONVERTEROBJS = radix51.o converter.o
MAIN = proj

.PHONY: depend clean

all:    test benchmark

# $(MAIN): $(OBJS)
# 	$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS)

test:	$(TESTOBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o test $(addprefix $(BUILDDIR)/, $(TESTOBJS)) $(LFLAGS)

converter:	$(CONVERTEROBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o converter $(addprefix $(BUILDDIR)/, $(CONVERTEROBJS))

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $(addprefix $(BUILDDIR)/, $@) $(LFLAGS)

radix51.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $(RADIXSRC) -o $(addprefix $(BUILDDIR)/, $@)

scalarmult.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c scalarmult.cpp -o $(addprefix $(BUILDDIR)/, $@)

tests.o: testsrc/tests.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $(BUILDDIR)/$@ $(LFLAGS)

converter.o: testsrc/converter.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $(BUILDDIR)/$@ $(LFLAGS)

clean:
	$(RM) $(BUILDDIR)/*.o *~ test converter

depend: $(BENCHSRCS)
	makedepend $(INCLUDES) $^
# DO NOT DELETE THIS LINE -- make depend needs it
