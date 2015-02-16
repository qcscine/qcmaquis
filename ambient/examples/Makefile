# Makefile example (icc)

SHELL := /bin/bash
CXXFLAGS?= -std=c++11 -O3 -DNDEBUG -DAMBIENT_GLOBALS -DAMBIENT_DISABLE_MPI -I../src -mkl

.DEFAULT: all
.PHONY: all clean

EXAMPLES = diffusion partitioned_vector vector atomic matrix mp_gemm_simple mp_gemm_log groups
all: $(EXAMPLES)

%: %.cpp
	$(CXX) -o $(basename $<) $< $(CXXFLAGS)

clean:
	-rm -f $(EXAMPLES)
