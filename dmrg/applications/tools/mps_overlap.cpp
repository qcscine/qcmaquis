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

#include "dmrg/sim/matrix_types.h"
#include "dmrg/block_matrix/indexing.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/contractions.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/mp_tensors/mpo_ops.h"

#if defined(USE_TWOU1)
typedef TwoU1 grp;
#elif defined(USE_TWOU1PG)
typedef TwoU1PG grp;
#elif defined(USE_SU2U1)
typedef SU2U1 grp;
#elif defined(USE_SU2U1PG)
typedef SU2U1PG grp;
#elif defined(USE_NONE)
typedef TrivialGroup grp;
#elif defined(USE_U1)
typedef U1 grp;
#endif


int main(int argc, char ** argv)
{
    try {
        if (argc != 3) {
            std::cout << "Usage: " << argv[0] << " <mps1.h5> <mps2.h5>" << std::endl;
            return 1;
        }
        MPS<matrix, grp> mps1, mps2;
        load(argv[1], mps1);
        load(argv[2], mps2);
        
        if (true) {
            std::cout << "<mps1 | mps2> = "                        << std::endl;
            std::cout                       << overlap(mps1, mps2) << std::endl;
            std::cout << "<mps2 | mps1> = "                        << std::endl;
            std::cout                       << overlap(mps2, mps1) << std::endl;
        }
        
        if (false) {
            operator_selector<matrix, grp>::type ident;
            for (int i=0; i<mps1.site_dim(0).size(); ++i)
                ident.insert_block(matrix::identity_matrix(mps1.site_dim(0)[i].second),
                                   mps1.site_dim(0)[i].first, mps1.site_dim(0)[i].first);
            
            MPO<matrix, grp> mpo;
            
            MPOTensor<matrix, grp> mpot;
            mpot.set(0,0, ident);
            mpo = MPO<matrix, grp>(mps1.length());
            for (int p=0; p<mps1.length(); ++p)
                mpo[p] = mpot;
            
            std::cout << "<mps1 | 1 | mps2> = " << expval(mps1, mps2, mpo) << std::endl;
            std::cout << "<mps2 | 1 | mps1> = " << expval(mps2, mps1, mpo) << std::endl;
        }
    } catch (std::exception& e) {
        std::cerr << "Error:" << std::endl << e.what() << std::endl;
        return 1;
    }
}
