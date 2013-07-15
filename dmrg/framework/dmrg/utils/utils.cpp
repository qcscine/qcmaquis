#include "utils.hpp"
#include "random.hpp"
#include "dmrg/utils/archive.h"
#include "dmrg/utils/logger.h"

// Init Comparison object
double cmp_with_prefactor::prefactor = 1.;

// Init random
dmrg_random::engine_t dmrg_random::engine = dmrg_random::engine_t(42);

// Init logger
namespace storage {
    Logger<storage::archive> log;
}
