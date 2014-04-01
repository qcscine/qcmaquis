# Makefile
SHELL := /bin/bash
BOOST_INCLUDE_PATH?=/project/h07/boost_1_53_0
AMBIENT_PATH?=/users/akosenko/maquis2014/src/ambient/src
AMBIENT_INCLUDE_PATH?=$(AMBIENT_PATH)
AMBIENT_LIB_PATH?=.
CXXFLAGS?= -std=c++11 -O3 -DNDEBUG

# icc:      export CXX=icpc        # -std=c++11 -O3 -DNDEBUG
# clang:    export CXX=g++         # -std=c++11 -O3 -DNDEBUG -DAMBIENT_SERIAL -Wno-array-bounds (-DAMBIENT_NOP_CHANNEL)
# gcc:      export CXX=g++-mp-4.7  # -std=c++11 -O3 -DNDEBUG -fopenmp (-DAMBIENT_NOP_CHANNEL)
# openmpi:  export OMPI_CXX=$CXX; export CXX=openmpicxx;
# mvapich:  export CXX=mpicxx;

# Notes:
# - for alps::hdf5 support include -I$ALPS_INCLUDE_PATH
# - add -lboost -L$BOOST_LIB_PATH -Wl,-rpath,$BOOST_LIB_PATH if need to link against boost
# - to change threading selection from auto - define CILK/OMP/SERIAL with AMBIENT_ prefix

.DEFAULT: all
.PHONY: all

all: libambient.a diffusion vector matrix

diffusion: libambient.a diffusion.cpp
	$(CXX) diffusion.cpp -o diffusion $(CXXFLAGS) -lambient -L$(AMBIENT_LIB_PATH) -I$(AMBIENT_INCLUDE_PATH)

vector: libambient.a vector.cpp
	$(CXX) vector.cpp -o vector $(CXXFLAGS) -lambient -L$(AMBIENT_LIB_PATH) -I$(AMBIENT_INCLUDE_PATH)

matrix: libambient.a matrix.cpp
	$(CXX) matrix.cpp  -o matrix  $(CXXFLAGS) -lambient -L$(AMBIENT_LIB_PATH) -I$(AMBIENT_INCLUDE_PATH) -mkl -I$(BOOST_INCLUDE_PATH)

libambient.a: $(AMBIENT_PATH)/ambient/ambient.cpp
	$(CXX) -c $(AMBIENT_PATH)/ambient/ambient.cpp $(CXXFLAGS) -I$(AMBIENT_PATH) -o ambient.cpp.o
	ar cr libambient.a ambient.cpp.o; ranlib libambient.a;
	rm -f ambient.cpp.o	

.PHONY: clean
clean:
	-rm -f libambient.a
	-rm -f diffusion
	-rm -f vector
	-rm -f matrix
