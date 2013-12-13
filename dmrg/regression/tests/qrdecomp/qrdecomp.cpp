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

#define BOOST_TEST_MAIN

#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <iterator>
#include <iostream>


using std::cerr;
using std::cout;
using std::endl;

#include "dmrg/block_matrix/detail/alps.hpp"

#include "dmrg/utils/DmrgParameters2.h"

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mps_initializers.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/models/coded/lattice.hpp"


typedef U1 SymmGroup;
typedef alps::numeric::matrix<double> matrix;

BOOST_AUTO_TEST_CASE( test )
{
    
    typedef std::vector<block_matrix<matrix, SymmGroup> > op_vec;
    int Nrep = 6;
    int M = 50;
    int L = 16;

    boost::shared_ptr<lattice_impl> lat_ptr(new ChainLattice(L));
    Lattice lattice(lat_ptr);
    
    DmrgParameters parms;
    parms.set("max_bond_dimension", M);

    // Bosons with Nmax=2
    Index<SymmGroup> phys;
    phys.insert(std::make_pair(0, 1));
    phys.insert(std::make_pair(1, 1));
    phys.insert(std::make_pair(2, 1));
    SymmGroup::charge initc = L/2;
    
    block_matrix<matrix, SymmGroup> ident = identity_matrix<matrix>(phys);
    op_vec ids_vec(1,ident);

    default_mps_init<matrix, SymmGroup> initializer(parms, std::vector<Index<SymmGroup> >(1, phys), initc, std::vector<int>(L,0));
    
    MPS<matrix,SymmGroup> mps;
    mps.resize(L); initializer(mps);
    double onorm = norm(mps);
    
    MPS<matrix, SymmGroup> mps1(mps), mps2(mps);
    
    // here canonize with SVD
    mps1.canonize(L/2, SVD);
    
    // here canonize with QR
    mps2.canonize(L/2, QR);
    
    { // measure norm
        double orig = norm(mps);
        double meas_mps1 = norm(mps1);
        double meas_mps2 = norm(mps2);

        maquis::cout << "Norm: ";
        maquis::cout << orig << " " << meas_mps1 << " " << meas_mps2 << std::endl;
    }
    
    { // measure density (at some positions)
        block_matrix<matrix, SymmGroup> op;
        op.insert_block(matrix(1,1,1), 1,1);
        op.insert_block(matrix(1,1,2), 2,2);
        
        for (int i=0; i<10; ++i) {
            int p = dmrg_random::uniform() * L;
            
            generate_mpo::MPOMaker<matrix, SymmGroup> mpom(lattice, ids_vec, ids_vec);
            generate_mpo::Operator_Term<matrix, SymmGroup> term;
            term.operators.push_back( std::make_pair(p, op) );
            term.fill_operator = ident;
            mpom.add_term(term);
            MPO<matrix, SymmGroup> mpo = mpom.create_mpo();
            
            double orig = expval(mps, mpo);
            double meas_mps1 = expval(mps1, mpo);
            double meas_mps2 = expval(mps2, mpo);
            
            maquis::cout << "Density at site " << p << ": ";
            maquis::cout << orig / onorm << " " << meas_mps1 << " " << meas_mps2 << std::endl;
        }
    }

    { // measure bond terms
        block_matrix<matrix, SymmGroup> op1, op2;
        op1.insert_block(matrix(1,1,1), 0,1);
        op1.insert_block(matrix(1,1,sqrt(2)), 1,2);
        op2.insert_block(matrix(1,1,1), 1,0);
        op2.insert_block(matrix(1,1,sqrt(2)), 2,1);
        
        generate_mpo::MPOMaker<matrix, SymmGroup> mpom(lattice, ids_vec, ids_vec);
        for (int p=0; p<L-1; ++p) {
            {
                generate_mpo::Operator_Term<matrix, SymmGroup> term;
                term.operators.push_back( std::make_pair(p, -1.*op1) );
                term.operators.push_back( std::make_pair(p+1, op2) );
                term.fill_operator = ident;
                mpom.add_term(term);
            }
            {
                generate_mpo::Operator_Term<matrix, SymmGroup> term;
                term.operators.push_back( std::make_pair(p, -1.*op2) );
                term.operators.push_back( std::make_pair(p+1, op1) );
                term.fill_operator = ident;
                mpom.add_term(term);
            }
        }
        MPO<matrix, SymmGroup> mpo = mpom.create_mpo();
        
        double orig = expval(mps, mpo);
        double meas_mps1 = expval(mps1, mpo);
        double meas_mps2 = expval(mps2, mpo);
        
        maquis::cout << "Bond energy: ";
        maquis::cout << orig / onorm << " " << meas_mps1 << " " << meas_mps2 << std::endl;
    }

    { // measure using canonization
        block_matrix<matrix, SymmGroup> op;
        op.insert_block(matrix(1,1,1), 1,1);
        op.insert_block(matrix(1,1,2), 2,2);
        
        for (int i=0; i<10; ++i) {
            int p = dmrg_random::uniform() * L;
            
            generate_mpo::MPOMaker<matrix, SymmGroup> mpom(lattice, ids_vec, ids_vec);
            generate_mpo::Operator_Term<matrix, SymmGroup> term;
            term.operators.push_back( std::make_pair(p, op) );
            term.fill_operator = ident;
            mpom.add_term(term);
            MPO<matrix, SymmGroup> mpo = mpom.create_mpo();
            
            double orig = expval(mps, mpo);
            
            // here canonize with SVD
            mps1.canonize(p, SVD);
            
            double meas_mps1 = mps1[p].scalar_overlap(contraction::local_op(mps1[p], op));
            
            // here canonize with QR
            mps2.canonize(p, QR);
            double meas_mps2 = mps2[p].scalar_overlap(contraction::local_op(mps2[p], op));
          
            double ref_value = orig / onorm;
            BOOST_CHECK_CLOSE(meas_mps1, meas_mps2, 0.00000001 );
            BOOST_CHECK_CLOSE(ref_value, meas_mps1, 0.00000001);
            BOOST_CHECK_CLOSE(ref_value, meas_mps2, 0.00000001 );
            
            maquis::cout << "Density at site " << p << ": ";
            maquis::cout << orig / onorm << " " << meas_mps1 << " " << meas_mps2 << std::endl;
        }
    }
    
}
