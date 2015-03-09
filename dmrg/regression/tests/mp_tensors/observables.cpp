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

#define BOOST_TEST_MAIN

#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <iterator>
#include <iostream>


using std::cerr;
using std::cout;
using std::endl;

#include "dmrg/block_matrix/detail/alps.hpp"

#include "dmrg/utils/DmrgParameters.h"

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mps_initializers.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/models/coded/lattice.hpp"


typedef U1 SymmGroup;
typedef alps::numeric::matrix<double> matrix;

double measure_local_expval(MPS<matrix, SymmGroup> const & mps, typename operator_selector<matrix, SymmGroup>::type const & ident,
                            typename operator_selector<matrix, SymmGroup>::type const & op, size_t pos)
{
    typedef typename operator_selector<matrix, SymmGroup>::type op_t;
    typedef std::vector<op_t> op_vec;
    boost::shared_ptr<lattice_impl> lat_ptr(new ChainLattice(mps.length()));
    Lattice lattice(lat_ptr);
    
    generate_mpo::MPOMaker<matrix, SymmGroup> mpom(lattice, op_vec(1,ident), op_vec(1,ident));
    generate_mpo::OperatorTerm<matrix, SymmGroup> term;
    term.operators.push_back( std::make_pair(pos, op) );
    term.fill_operator = ident;
    mpom.add_term(term);
    MPO<matrix, SymmGroup> mpo = mpom.create_mpo();
    
    return expval(mps, mpo);
}

double measure_local_overlap(MPS<matrix, SymmGroup> const & mps,
                             typename operator_selector<matrix, SymmGroup>::type const & op, size_t pos)
{
    // asuming mps is canonized at site pos!
    return mps[pos].scalar_overlap(contraction::local_op(mps[pos], op));
}

double measure_local_trace(MPS<matrix, SymmGroup> const & mps,
                           block_matrix<matrix, SymmGroup> const & op, size_t pos)
{
    // asuming mps is canonized at site pos!
    block_matrix<matrix, SymmGroup> dm = contraction::density_matrix(mps[pos], mps[pos]);
    block_matrix<matrix, SymmGroup> tmp;
    gemm(op, dm, tmp);
    return trace(tmp);
}

double measure_local_trace_2(MPS<matrix, SymmGroup> const & mps,
                             block_matrix<matrix, SymmGroup> const & op, size_t pos)
{
    // asuming mps is canonized at site pos!
    block_matrix<matrix, SymmGroup> dm = contraction::density_matrix_2(mps[pos], mps[pos]);
    block_matrix<matrix, SymmGroup> tmp;
    gemm(op, dm, tmp);
    return trace(tmp);
}

BOOST_AUTO_TEST_CASE( obs_bosons_nmax2 )
{
    typedef typename operator_selector<matrix, SymmGroup>::type op_t;
    typedef block_matrix<matrix, SymmGroup> bm_t;
    
    int Nrep = 8;
    int M = 20;
    int L = 16;

    DmrgParameters parms;
    parms.set("max_bond_dimension", M);
    
    // Bosons with Nmax=2
    Index<SymmGroup> phys;
    phys.insert(std::make_pair(0, 1));
    phys.insert(std::make_pair(1, 1));
    phys.insert(std::make_pair(2, 1));
    SymmGroup::charge initc = L/2;
    
    op_t densop;
    densop.insert_block(matrix(1,1,1), 1,1);
    densop.insert_block(matrix(1,1,2), 2,2);
    op_t ident = identity_matrix<op_t>(phys);

    bm_t densop_bm;
    densop_bm.insert_block(matrix(1,1,1), 1,1);
    densop_bm.insert_block(matrix(1,1,2), 2,2);
    
    default_mps_init<matrix, SymmGroup> initializer(parms, std::vector<Index<SymmGroup> >(1, phys), initc, std::vector<int>(L,0));
    
    MPS<matrix,SymmGroup> mps;
    mps.resize(L); initializer(mps);
    mps.normalize_left();
    
    { // measure density (at some positions)
        
        for (int i=0; i<Nrep; ++i) {
            int p = dmrg_random::uniform() * L;
            mps.canonize(p);
            
            double meas[4];
            meas[0] = measure_local_expval(mps, ident, densop, p);
            meas[1] = measure_local_overlap(mps, densop, p);
            meas[2] = measure_local_trace(mps, densop_bm, p);
            meas[3] = measure_local_trace_2(mps, densop_bm, p);
            
            for (int i=1; i<4; ++i)
                BOOST_CHECK_CLOSE(meas[0], meas[i], 1e-8 );
        }
        
    }
    
}
