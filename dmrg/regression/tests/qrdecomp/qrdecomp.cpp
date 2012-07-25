
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


int main() {
    
    int Nrep = 10;
    int M = 50;
    int L = 40;
    
    // Bosons with Nmax=2
    Index<SymmGroup> phys;
    phys.insert(std::make_pair(0, 1));
    phys.insert(std::make_pair(1, 1));
    phys.insert(std::make_pair(2, 1));
    SymmGroup::charge initc = L/2;
    
    block_matrix<matrix, SymmGroup> ident = identity_matrix<matrix>(phys);
    
    default_mps_init<matrix, SymmGroup> initializer;
    
    MPS<matrix,SymmGroup> mps;
    mps.resize(L); initializer(mps, M, phys, initc);
    double onorm = norm(mps);
    
    MPS<matrix, SymmGroup> mps1(mps), mps2(mps);
    
    // todo: here set default decomposition to SVD
    mps1.canonize(20);
    
    // todo: here set default decomposition to QR
    mps2.canonize(20);
    
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
            
            generate_mpo::MPOMaker<matrix, SymmGroup> mpom(L, ident);
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
        
        generate_mpo::MPOMaker<matrix, SymmGroup> mpom(L, ident);
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
            
            generate_mpo::MPOMaker<matrix, SymmGroup> mpom(L, ident);
            generate_mpo::Operator_Term<matrix, SymmGroup> term;
            term.operators.push_back( std::make_pair(p, op) );
            term.fill_operator = ident;
            mpom.add_term(term);
            MPO<matrix, SymmGroup> mpo = mpom.create_mpo();
            
            double orig = expval(mps, mpo);
            
            // todo: here set default decomposition to SVD
            mps1.canonize(p);
            
            double meas_mps1 = mps1[p].scalar_overlap(contraction::local_op(mps1[p], op));
            
            // todo: here set default decomposition to QR
            mps2.canonize(p);
            double meas_mps2 = mps2[p].scalar_overlap(contraction::local_op(mps2[p], op));
            
            maquis::cout << "Density at site " << p << ": ";
            maquis::cout << orig / onorm << " " << meas_mps1 << " " << meas_mps2 << std::endl;
        }
    }
    
}
