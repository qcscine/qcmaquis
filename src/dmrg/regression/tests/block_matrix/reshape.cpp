#define BOOST_TEST_MODULE example

#include <cmath>
#include <iterator>
#include <iostream>

#include <boost/test/included/unit_test.hpp>
#include <boost/geometry/geometries/adapted/boost_array.hpp>

using std::cerr;
using std::cout;
using std::endl;

//#include <boost/filesystem.hpp>

#include "alps/numeric/matrix.hpp"
#include "alps/numeric/matrix/algorithms.hpp"
#include "dmrg/kernels/alps_matrix.hpp"
#include "alps/numeric/diagonal_matrix.hpp"

//#include <alps/hdf5.hpp>

#include "dmrg/utils/DmrgParameters2.h"

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/contractions.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/mp_tensors/mpo_ops.h"
#include "dmrg/mp_tensors/mps_initializers.h"

#include "dmrg/models/generate_mpo.hpp"

//#include "dmrg/utils/logger.h"
//#include "dmrg/utils/random.hpp"
//#include "utils/timings.h"

typedef U1 SymmGroup;
typedef alps::numeric::matrix<double> Matrix;


std::vector<double> density(MPS<Matrix, SymmGroup> const & mps, block_matrix<Matrix, SymmGroup> const & dens_op,
                            Index<SymmGroup> const & phys, block_matrix<Matrix, SymmGroup> const & ident)
{
    std::vector<double> vals(mps.length());
    
    for (int p=0; p<mps.length(); p++) {
        generate_mpo::MPOMaker<Matrix, SymmGroup> mpom(mps.length(), ident);
        generate_mpo::Operator_Term<Matrix, SymmGroup> term;
        term.operators.push_back( std::make_pair(p, dens_op) );
        term.fill_operator = ident;
        mpom.add_term(term);
        MPO<Matrix, SymmGroup> mpo = mpom.create_mpo();
        
        double val = expval(mps, mpo);
        vals[p] = val;
    }
    return vals;
}

MPS<Matrix, SymmGroup> state_mps(std::vector<short> const & state, Index<SymmGroup> const & phys)
{
    MPS<Matrix, SymmGroup> mps(state.size());
    
    Index<SymmGroup> curr_i;
    curr_i.insert(std::make_pair(0, 1));
    for (int i=0; i<state.size(); ++i)
    {
        Index<SymmGroup> new_i;
        new_i.insert(std::make_pair(curr_i[0].first + state[i], 1));
        mps[i] = MPSTensor<Matrix, SymmGroup>(phys, curr_i, new_i, false, 1);
        curr_i = new_i;
    }
    return mps;
}

std::vector<std::vector<short> > create_basis (int L, int N)
{
    std::vector<std::vector<short> > basis;
    for (int n2=0; n2*2<=N; ++n2)
    {
        int n1 = N - 2*n2;
        int n0 = std::max(0, L - n1 - n2);
        std::vector<short> state(L, 0);
        for (int i=n0; i<n0+n1; ++i)
            state[i] = 1;
        for (int i=n0+n1; i<n0+n1+n2; ++i)
            state[i] = 2;
        do {
            basis.push_back(state);
        } while ( std::next_permutation(state.begin(),state.end()) );
    }
    
    return basis;
}

BOOST_AUTO_TEST_CASE( free_test_function ) {
  // it is just a test,horrible 
/*
    std::vector<double> res;  
    std::vector<double> ref = { 0.224538,     
                         0.184199,     
                         0.227943,     
                         0.159481,    
                         0.205156,     
                         0.143854,     
                         0.184408,     
                         0.196194,     
                         0.102792,     
                         0.103562,     
                         0.148949,     
                         0.140995,     
                         0.0607287,     
                         0.181595 ,    
                         0.250387 ,    
                         0.497258 ,   
                         0.0563364 ,    
                         0.194877  ,  
                         0.687507  ,  
                         1.0492
                        }; 
 
    int Nrep = 1;
    
    int N = 4;
    int L = 5*N;
    SymmGroup::charge initc = L/N;
    int M = 600;
    
    // Bosons with Nmax=2
    Index<SymmGroup> phys;
    phys.insert(std::make_pair(0, 1));
    phys.insert(std::make_pair(1, 1));
    phys.insert(std::make_pair(2, 1));
    
    // Identity operator
    block_matrix<Matrix, SymmGroup> ident;
    ident.insert_block(Matrix(1,1, 1), 0, 0);
    ident.insert_block(Matrix(1,1, 1), 1, 1);
    ident.insert_block(Matrix(1,1, 1), 2, 2);
    
    // Density operator
    block_matrix<Matrix, SymmGroup> op;
    op.insert_block(Matrix(1,1, 0), 0, 0);
    op.insert_block(Matrix(1,1, 1), 1, 1);
    op.insert_block(Matrix(1,1, 2), 2, 2);
    
    //    const_mps_init<Matrix, SymmGroup> initializer;
    //    const_beta_mps_init<Matrix, SymmGroup> initializer;
    //    default_mps_init<Matrix, SymmGroup> initializer;
    //    linear_const_mps_init<Matrix> initializer;
        linear_mps_init<Matrix> initializer;
    //linear_beta_const_mps_init<Matrix> initializer;
    
    
    //    std::vector<std::vector<short> > basis = create_basis(L, initc);
    
    std::vector<std::vector<double> > dens(Nrep, std::vector<double>(L));
    //    std::vector<std::vector<double> > over(Nrep, std::vector<double>(basis.size()));
    for (int i=0; i<Nrep; i++) {
        MPS<Matrix, SymmGroup> mps(L, M, phys, initc, initializer);
        if (false) {
            for (int k=0; k<L; ++k) {
                std::cout << "* MPS site " << k << ":" << std::endl;
                std::cout << mps[k];
            }
        }
        dens[i] = density(mps, op, phys, ident);
        //        for (int j=0; j<basis.size(); ++j)
        //            over[i][j] = overlap(state_mps(basis[j], phys), mps);
        std::cout << "# Measurement " << i << " done." << std::endl;
    }
    
    //    std::cout << "Length of basis: " << basis.size() << std::endl;
    std::cout << "== DENSITIES ==" << std::endl;
    for (int p=0; p<L; p++)
    {
        for (int i=0; i<Nrep; i++)
            res[i] = dens[i][p];
    }
    
    //    if (false) {
    //        std::cout << "== OVERLAPS ==" << std::endl;
    //        for (int j=0; j<basis.size(); j++)
    //        {
    //            printf("#%03d: | ", j+1);
    //            std::copy(basis[j].begin(), basis[j].end(), std::ostream_iterator<short>(std::cout, " "));
    //            std::cout << ">\t";
    //            for (int i=0; i<Nrep; i++)
    //                std::cout << over[i][j] << "\t";
    //            std::cout << std::endl;
    //        }
    //    }
   BOOST_CHECK_EQUAL(res,ref);
*/
}
