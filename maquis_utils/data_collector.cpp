
#include "utils/data_collector.hpp"

DCOLLECTOR_CREATE(gemm_collector, "GEMM sizes", 100)
DCOLLECTOR_CREATE(svd_collector, "SVD sizes", 100)
DCOLLECTOR_CREATE(num_blocks_gemm_collector, "GEMM num blocks", 50)
DCOLLECTOR_CREATE(num_blocks_svd_collector, "SVD num blocks", 50)

