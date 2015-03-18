/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
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

#include <iostream>
#include <fstream>

#include "dmrg/utils/DmrgParameters.h"

#include "dmrg/sim/symm_factory.h"

int main(int argc, char ** argv)
{
    try {
        if (argc != 3)
        {
            maquis::cout << "Usage: <parms> <model_parms>" << std::endl;
            exit(1);
        }
        
        maquis::cout.precision(10);
        
        /// Loading parameters
        std::ifstream param_file(argv[1]);
        if (!param_file) {
            maquis::cerr << "Could not open parameter file." << std::endl;
            exit(1);
        }
        DmrgParameters parms(param_file);
        
        /// Loading model
        std::ifstream model_file(argv[2]);
        if (!model_file) {
            maquis::cerr << "Could not open model file." << std::endl;
            exit(1);
        }
        ModelParameters model_parms(model_file);
        
        
        /// Invoke symmetry factory
        maquis::dmrg::symm_factory(parms, model_parms);
        
        
    } catch (std::exception & e) {
        maquis::cerr << "Exception thrown!" << std::endl;
        maquis::cerr << e.what() << std::endl;
        exit(1);
    }
}

