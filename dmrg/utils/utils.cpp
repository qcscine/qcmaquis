#ifdef MPI_PARALLEL
#include "mpi.h"
#endif

#include "utils.hpp"
#include "zout.hpp"
#include "random.hpp"

// Init Comparison object
double cmp_with_prefactor::prefactor = 1.;

// Init zout
master_cout zout;
group_master_cout gzout;

// Init random
dmrg_random::engine_t dmrg_random::engine = dmrg_random::engine_t(42);

