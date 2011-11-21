#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

#include "ambient/utils/zout.hpp"

namespace ambient {
    // Init zout
    master_cout zout;
    group_master_cout gzout;
}
