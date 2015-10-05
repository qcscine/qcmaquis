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

#ifdef USE_AMBIENT
#include <mpi.h>
#endif
#include <cmath>
#include <iterator>
#include <iostream>
#include <sys/time.h>
#include <sys/stat.h>

using std::cerr;
using std::cout;
using std::endl;

#ifdef USE_AMBIENT
    #include "dmrg/block_matrix/detail/ambient.hpp"
    typedef ambient::tiles<ambient::matrix<double> > matrix;
#else
    #include "dmrg/block_matrix/detail/alps.hpp"
    typedef alps::numeric::matrix<double> matrix;
#endif

#include "dmrg/block_matrix/block_matrix.h"

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/mp_tensors/mps.h"

#if defined(USE_NONE)
typedef TrivialGroup grp;
#elif defined(USE_U1)
typedef U1 grp;
#elif defined(USE_TWOU1)
typedef TwoU1 grp;
#else
#error "SymmGroup has to be defined explicitly. (-DUSE_NONE, -DUSE_U1, -DUSE_TWOU1)"
#endif


int main(int argc, char ** argv)
{
    try {
        if (argc != 3) {
            std::cout << "Usage: " << argv[0] << " <in.h5> <out.chkp>" << std::endl;
            return 1;
        }
        
        alps::hdf5::archive ih5(argv[1]);
        
        /// hard-coded old format
        std::vector<std::string> children = ih5.list_children("/state/MPS");
        MPS<matrix, grp> mps(children.size());
        for (std::vector<std::string>::const_iterator it = children.begin(); it != children.end(); ++it)
            ih5["/state/MPS/" + *it] >> mps[alps::cast<std::size_t>(*it)];
        
        /// create chkp directory
        if (!boost::filesystem::exists(argv[2]))
            boost::filesystem::create_directory(argv[2]);
        
        /// save in new format
        save(argv[2], mps);
        
        /// copy /version, /parameters and /status
        alps::hdf5::archive oh5(argv[2] + std::string("/props.h5"), "w");
        
        std::map<std::string, int> status;
        ih5["/status"] >> status;
        oh5["/status"] << status;

        std::string version;
        ih5["/version"] >> version;
        oh5["/version"] << version;
        
//        alps::ngs::params parameters;
//        ih5["/parameters"] >> parameters;
//        oh5["/parameters"] << parameters;

        
    } catch (std::exception& e) {
        std::cerr << "Error:" << std::endl << e.what() << std::endl;
        return 1;
    }
}
