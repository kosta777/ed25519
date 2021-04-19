mkdir -p StdC-optimize/build-O0
mkdir -p StdC-optimize/build-ntv
mkdir -p StdC-optimize/build-O3
mkdir -p SSE/build-O0
mkdir -p SSE/build-ntv
mkdir -p SSE/build-O3

make clean BUILDSUBDIR=build-O0     BINDIR=StdC-optimize/build-O0/
make clean BUILDSUBDIR=build-ntv    BINDIR=StdC-optimize/build-ntv/
make clean BUILDSUBDIR=build-O3     BINDIR=StdC-optimize/build-O3/
make clean BUILDSUBDIR=build-O0     BINDIR=SSE/build-O0/
make clean BUILDSUBDIR=build-ntv    BINDIR=SSE/build-ntv/
make clean BUILDSUBDIR=build-O3     BINDIR=SSE/build-O3/

make OPTFLAGS="-O0 -march=native"                       BUILDSUBDIR=build-O0    BINDIR=StdC-optimize/build-O0/  OPTDIR=StdC-optimize
make OPTFLAGS="-O3 -fno-tree-vectorize -march=native"   BUILDSUBDIR=build-ntv   BINDIR=StdC-optimize/build-ntv/ OPTDIR=StdC-optimize
make OPTFLAGS="-O3 -march=native"                       BUILDSUBDIR=build-O3    BINDIR=StdC-optimize/build-O3/  OPTDIR=StdC-optimize
make OPTFLAGS="-O0 -march=native"                       BUILDSUBDIR=build-O0    BINDIR=SSE/build-O0/            OPTDIR=SSE
make OPTFLAGS="-O3 -fno-tree-vectorize -march=native"   BUILDSUBDIR=build-ntv   BINDIR=SSE/build-ntv/           OPTDIR=SSE
make OPTFLAGS="-O3 -march=native"                       BUILDSUBDIR=build-O3    BINDIR=SSE/build-O3/            OPTDIR=SSE