#ifdef MPI_PARALLEL
#include "mpi.h"
#endif

#include "utils.hpp"
#include "random.hpp"

// Init Comparison object
double cmp_with_prefactor::prefactor = 1.;

// Init random
dmrg_random::engine_t dmrg_random::engine = dmrg_random::engine_t(42);

