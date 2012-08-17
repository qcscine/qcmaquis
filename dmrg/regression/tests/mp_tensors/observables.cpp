#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
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


typedef U1 SymmGroup;
typedef alps::numeric::matrix<double> matrix;

double measure_local_expval(MPS<matrix, SymmGroup> const & mps, block_matrix<matrix, SymmGroup> const & ident,
                            block_matrix<matrix, SymmGroup> const & op, size_t pos)
{
    generate_mpo::MPOMaker<matrix, SymmGroup> mpom(mps.length(), ident);
    generate_mpo::Operator_Term<matrix, SymmGroup> term;
    term.operators.push_back( std::make_pair(pos, op) );
    term.fill_operator = ident;
    mpom.add_term(term);
    MPO<matrix, SymmGroup> mpo = mpom.create_mpo();
    
    return expval(mps, mpo);
}

double measure_local_overlap(MPS<matrix, SymmGroup> const & mps,
                             block_matrix<matrix, SymmGroup> const & op, size_t pos)
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
    
    int Nrep = 10;
    int M = 50;
    int L = 40;
    
    // Bosons with Nmax=2
    Index<SymmGroup> phys;
    phys.insert(std::make_pair(0, 1));
    phys.insert(std::make_pair(1, 1));
    phys.insert(std::make_pair(2, 1));
    SymmGroup::charge initc = L/2;
    
    block_matrix<matrix, SymmGroup> densop;
    densop.insert_block(matrix(1,1,1), 1,1);
    densop.insert_block(matrix(1,1,2), 2,2);
    
    block_matrix<matrix, SymmGroup> ident = identity_matrix(block_matrix<matrix,SymmGroup>(),phys);
    
    default_mps_init<matrix, SymmGroup> initializer;
    
    MPS<matrix,SymmGroup> mps;
    mps.resize(L); initializer(mps, M, phys, initc);
    mps.normalize_left();
    
    { // measure density (at some positions)
        
        for (int i=0; i<Nrep; ++i) {
            int p = dmrg_random::uniform() * L;
            mps.canonize(p);
            
            double meas[4];
            meas[0] = measure_local_expval(mps, ident, densop, p);
            meas[1] = measure_local_overlap(mps, densop, p);
            meas[2] = measure_local_trace(mps, densop, p);
            meas[3] = measure_local_trace_2(mps, densop, p);
            
            for (int i=1; i<4; ++i)
                BOOST_CHECK_CLOSE(meas[0], meas[i], 1e-8 );
        }
        
    }
    
}
