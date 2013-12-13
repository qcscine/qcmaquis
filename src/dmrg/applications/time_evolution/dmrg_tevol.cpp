/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 * 
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 * 
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not, the license is also
 * available from http://alps.comp-phys.org/.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
 * DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

#include "dmrg/sim/symm_factory.h"

#include <cmath>
#include <iterator>
#include <iostream>
#include <sys/stat.h>

#include "utils/data_collector.hpp"
#include "dmrg/utils/DmrgParameters2.h"
#include "utils/timings.h"

#include <boost/function.hpp>
#undef tolower
#undef toupper
#include <boost/tokenizer.hpp>

int main(int argc, char ** argv)
{
    if (argc != 2 && argc != 3)
    {
        maquis::cout << "Usage: <parms> [<model_parms>]" << std::endl;
        exit(1);
    }
    
    maquis::cout.precision(10);
    
    /// Load parameters
    std::ifstream param_file(argv[1]);
    if (!param_file)
        throw std::runtime_error("Could not open parameter file.");
    DmrgParameters parms(param_file);
    
    /// Load model parameters from second input (if needed)
    std::string model_file;
    if (parms.is_set("model_file")) model_file = parms["model_file"].str();
    if (argc == 3)                  model_file = std::string(argv[2]);
    if (!model_file.empty()) {
        std::ifstream model_ifs(model_file.c_str());
        if (!model_ifs)
            throw std::runtime_error("Could not open model_parms file.");
        parms << ModelParameters(model_ifs);
    }
    
    DCOLLECTOR_SET_SIZE(gemm_collector, parms["max_bond_dimension"]+1)
    DCOLLECTOR_SET_SIZE(svd_collector, parms["max_bond_dimension"]+1)
    
    timeval now, then, snow, sthen;
    gettimeofday(&now, NULL);
    
    try {
        maquis::dmrg::symm_factory(parms);
    } catch (std::exception & e) {
        maquis::cerr << "Exception thrown!" << std::endl;
        maquis::cerr << e.what() << std::endl;
        exit(1);
    }
    
    gettimeofday(&then, NULL);
    double elapsed = then.tv_sec-now.tv_sec + 1e-6 * (then.tv_usec-now.tv_usec);
    
    DCOLLECTOR_SAVE_TO_FILE(gemm_collector, "collectors.h5", "/results")
    DCOLLECTOR_SAVE_TO_FILE(svd_collector, "collectors.h5", "/results")
    
    maquis::cout << "Task took " << elapsed << " seconds." << std::endl;
}

