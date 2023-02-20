/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#include <cmath>
#include <iterator>
#include <iostream>
#include <sys/stat.h>

#include "utils/data_collector.hpp"

#include "dmrg/utils/DmrgOptions.h"
#include "dmrg/utils/DmrgParameters.h"
#include "utils/timings.h"

#include <boost/function.hpp>

#include "mg_dmrg_traits.h"
#include "dmrg/sim/symmetry_factory.h"

int main(int argc, char ** argv)
{
    DmrgOptions opt(argc, argv);
    if (opt.valid) {
        
        maquis::cout.precision(10);
        
        //    DCOLLECTOR_SET_SIZE(gemm_collector, parms["max_bond_dimension"]+1)
        //    DCOLLECTOR_SET_SIZE(svd_collector, parms["max_bond_dimension"]+1)
        
        timeval now, then, snow, sthen;
        gettimeofday(&now, NULL);
        
//        try{
        simulation_traits::shared_ptr sim = dmrg::symmetry_factory<simulation_traits>(opt.parms);
        sim->run(opt.parms);

//        }catch(std::exception & e){
//            maquis::cerr << "Exception thrown!" << std::endl;
//            maquis::cerr << e.what() << std::endl;
//            exit(1);
//        }
        
        gettimeofday(&then, NULL);
        double elapsed = then.tv_sec-now.tv_sec + 1e-6 * (then.tv_usec-now.tv_usec);
        
        //    DCOLLECTOR_SAVE_TO_FILE(gemm_collector, "collectors.h5", "/results")
        //    DCOLLECTOR_SAVE_TO_FILE(svd_collector, "collectors.h5", "/results")
        
        maquis::cout << "Task took " << elapsed << " seconds." << std::endl;
    }
}

