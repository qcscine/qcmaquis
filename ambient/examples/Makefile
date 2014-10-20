# Makefile
SHELL := /bin/bash
BOOST_INCLUDE_PATH?=/project/h07/boost_1_53_0
AMBIENT_PATH?=/users/akosenko/maquis2014/src/ambient/src
AMBIENT_INCLUDE_PATH?=$(AMBIENT_PATH)
AMBIENT_LIB_PATH?=.
COMMON?= -std=c++11 -O3 -DNDEBUG
CXXFLAGS?= $(COMMON)

# CHANNEL?= -DDISABLE_MPI # no mpi
# CXXFLAGS?= $(COMMON) $(CHANNEL) -Wfatal-errors -gcc-name=/opt/local/bin/gcc-mp-4.7 -gxx-name=/opt/local/bin/g++-mp-4.7 -gcc-version=473 -I/opt/local/include/gcc47/c++ -I/opt/local/include/gcc47/c++/x86_64-apple-darwin13/ -L/opt/local/lib/gcc47/        # icc
# CXXFLAGS?= $(COMMON) $(CHANNEL) -Wno-array-bounds -DAMBIENT_SERIAL # clang
# CXXFLAGS?= $(COMMON) $(CHANNEL) -fopenmp                           # gcc

# Notes:
# - for alps::hdf5 support include -I$ALPS_INCLUDE_PATH
# - add -lboost -L$BOOST_LIB_PATH -Wl,-rpath,$BOOST_LIB_PATH if need to link against boost
# - to change threading selection from auto - define CILK/OMP/SERIAL with AMBIENT_ prefix

.DEFAULT: all
.PHONY: all

all: libambient.a diffusion partitioned_vector vector atomic matrix mp_gemm_simple mp_gemm_log groups

diffusion: libambient.a diffusion.cpp
	$(CXX) diffusion.cpp -o diffusion $(CXXFLAGS) -lambient -L$(AMBIENT_LIB_PATH) -I$(AMBIENT_INCLUDE_PATH)

partitioned_vector: libambient.a partitioned_vector.cpp
	$(CXX) partitioned_vector.cpp -o partitioned_vector $(CXXFLAGS) -lambient -L$(AMBIENT_LIB_PATH) -I$(AMBIENT_INCLUDE_PATH)

vector: libambient.a vector.cpp
	$(CXX) vector.cpp -o vector $(CXXFLAGS) -lambient -L$(AMBIENT_LIB_PATH) -I$(AMBIENT_INCLUDE_PATH)

atomic: libambient.a atomic.cpp
	$(CXX) atomic.cpp -o atomic $(CXXFLAGS) -lambient -L$(AMBIENT_LIB_PATH) -I$(AMBIENT_INCLUDE_PATH)

matrix: libambient.a matrix.cpp
	$(CXX) matrix.cpp -o matrix $(CXXFLAGS) -lambient -L$(AMBIENT_LIB_PATH) -I$(AMBIENT_INCLUDE_PATH) -mkl -I$(BOOST_INCLUDE_PATH)

mp_gemm_simple: libambient.a mp_gemm_simple.cpp
	$(CXX) mp_gemm_simple.cpp -o mp_gemm_simple $(CXXFLAGS) -lambient -L$(AMBIENT_LIB_PATH) -I$(AMBIENT_INCLUDE_PATH) -mkl -I$(BOOST_INCLUDE_PATH)

mp_gemm_log: libambient.a mp_gemm_log.cpp
	$(CXX) mp_gemm_log.cpp -o mp_gemm_log $(CXXFLAGS) -lambient -L$(AMBIENT_LIB_PATH) -I$(AMBIENT_INCLUDE_PATH) -mkl -I$(BOOST_INCLUDE_PATH)

groups: libambient.a groups.cpp
	$(CXX) groups.cpp -o groups $(CXXFLAGS) -lambient -L$(AMBIENT_LIB_PATH) -I$(AMBIENT_INCLUDE_PATH) -mkl -I$(BOOST_INCLUDE_PATH)

libambient.a: $(AMBIENT_PATH)/ambient/ambient.cpp
	$(CXX) -c $(AMBIENT_PATH)/ambient/ambient.cpp $(CXXFLAGS) -I$(AMBIENT_PATH) -o ambient.cpp.o
	ar cr libambient.a ambient.cpp.o; ranlib libambient.a;
	rm -f ambient.cpp.o	

.PHONY: clean
clean:
	-rm -f libambient.a
	-rm -f diffusion
	-rm -f partitioned_vector
	-rm -f vector
	-rm -f atomic
	-rm -f matrix
	-rm -f mp_gemm_simple
	-rm -f mp_gemm_log
	-rm -f groups
