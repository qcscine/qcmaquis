# Makefile example (icc)

SHELL := /bin/bash
CXXFLAGS?= -std=c++11 -O3 -DNDEBUG -I../../include -DAMBIENT_MPI=MPI_DISABLE -mkl

.DEFAULT: all
.PHONY: all clean

EXAMPLES = diffusion partitioned_vector vector atomic matrix mp_gemm_simple mp_gemm_log groups future
all: $(EXAMPLES)

%: %.cpp
	$(CXX) -o $(basename $<) $< $(CXXFLAGS)

clean:
	-rm -f $(EXAMPLES)
