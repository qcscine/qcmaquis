#include <cmath>
#include <iterator>
#include <iostream>

using std::cerr;
using std::cout;
using std::endl;

#include "dmrg/block_matrix/detail/alps.hpp"

#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/super_mpo.h"


typedef alps::numeric::matrix<double> matrix;


template <class SymmGroup>
void mpo1_test(Index<SymmGroup> const& phys,
               block_matrix<matrix, SymmGroup> const& op1,
               block_matrix<matrix, SymmGroup> const& op2,
               block_matrix<matrix, SymmGroup> const& op3)
{
    maquis::cout << "TEST MPO 1" << std::endl;
    
    /// MPO with three sites.
    ///   - op1 - op1 - op2 - op3 -
    ///   `      `op2 - op3 - op1´
    ///    `op2 - op3 - op1´
    MPO<matrix,SymmGroup> mpo(4);
    int site = 0;
    {
        MPOTensor<matrix, SymmGroup> mpot(1, 2);
        mpot(0,0) = op1;
        mpot(0,1) = op2;
        
        mpo[site++] = mpot;
    }
    
    {
        MPOTensor<matrix, SymmGroup> mpot(2, 3);
        mpot(0,0) = op1;
        mpot(1,1) = op3;
        mpot(0,2) = op2;
        
        mpo[site++] = mpot;
    }
    
    {
        MPOTensor<matrix, SymmGroup> mpot(3, 3);
        mpot(1,1) = op1;
        mpot(2,1) = op3;
        mpot(0,2) = op2;
        
        mpo[site++] = mpot;
    }
    
    {
        MPOTensor<matrix, SymmGroup> mpot(3, 2);
        mpot(1,1) = op1;
        mpot(2,1) = op3;
        
        mpo[site++] = mpot;
    }
    
    
    MPS<matrix, SymmGroup> super_mps = mpo_to_smps(mpo, phys);
    std::cout << "MPS description:" << std::endl << super_mps.description();
}

template <class SymmGroup>
void mpo2_test(Index<SymmGroup> const& phys,
               block_matrix<matrix, SymmGroup> const& op1,
               block_matrix<matrix, SymmGroup> const& op2,
               block_matrix<matrix, SymmGroup> const& op3)
{
    maquis::cout << "TEST MPO 2" << std::endl;
    
    /// MPO with three sites.
    ///   - op1 - op1 - op2 -
    ///   `      `op2 - op3 -
    ///    `op2 - op3 - op1´
    MPO<matrix,SymmGroup> mpo(3);
    int site = 0;
    {
        MPOTensor<matrix, SymmGroup> mpot(1, 2);
        mpot(0,0) = op1;
        mpot(0,1) = op2;
        
        mpo[site++] = mpot;
    }
    
    {
        MPOTensor<matrix, SymmGroup> mpot(2, 3);
        mpot(0,0) = op1;
        mpot(1,1) = op3;
        mpot(0,2) = op2;
        
        mpo[site++] = mpot;
    }
    
    {
        MPOTensor<matrix, SymmGroup> mpot(3, 3);
        mpot(1,1) = op1;
        mpot(2,1) = op3;
        mpot(0,2) = op2;
        
        mpo[site++] = mpot;
    }
    
    
    MPS<matrix, SymmGroup> super_mps = mpo_to_smps(mpo, phys);
    std::cout << "MPS description:" << std::endl << super_mps.description();
}

template <class SymmGroup>
void mpo3_test(Index<SymmGroup> const& phys,
               block_matrix<matrix, SymmGroup> const& op1,
               block_matrix<matrix, SymmGroup> const& op2,
               block_matrix<matrix, SymmGroup> const& op3)
{
    maquis::cout << "TEST MPO 3" << std::endl;
    
    /// MPO with four sites.
    ///  - op1 - op3 - op1 - op1
    /// `            `           `
    ///  `            `op2 - op3 -`
    ///   `op3 - op1 ´            ´
    ///              ` op3 - op2 ´

    MPO<matrix,SymmGroup> mpo(4);
    int site = 0;
    {
        MPOTensor<matrix, SymmGroup> mpot(1, 2);
        mpot(0,0) = op1;
        mpot(0,1) = op3;
        
        mpo[site++] = mpot;
    }
    
    {
        MPOTensor<matrix, SymmGroup> mpot(2, 2);
        mpot(0,0) = op3;
        mpot(1,1) = op1;
        
        mpo[site++] = mpot;
    }
    
    {
        MPOTensor<matrix, SymmGroup> mpot(2, 3);
        // this mpo seems impossible...
        mpot(0,0) = op1;
        mpot(0,1) = op2;
        mpot(0,2) = op3;
        
        mpo[site++] = mpot;
    }

    {
        MPOTensor<matrix, SymmGroup> mpot(3, 3);
        mpot(0,1) = op1;
        mpot(1,1) = op3;
        mpot(2,1) = op2;
        
        mpo[site++] = mpot;
    }

    MPS<matrix, SymmGroup> super_mps = mpo_to_smps(mpo, phys);
    std::cout << "MPS description:" << std::endl << super_mps.description();
}


template <class SymmGroup>
void all_tests(Index<SymmGroup> const& phys,
               block_matrix<matrix, SymmGroup> const& op1,
               block_matrix<matrix, SymmGroup> const& op2,
               block_matrix<matrix, SymmGroup> const& op3)
{
    
    mpo1_test(phys, op1, op2, op3);
    mpo2_test(phys, op1, op2, op3);
}

void test_none()
{
    typedef TrivialGroup SymmGroup;
    maquis::cout << "TESTING NONE SYMMETRY" << std::endl;
    
    SymmGroup::charge C = SymmGroup::IdentityCharge;
    
    Index<SymmGroup> phys;
    phys.insert(std::make_pair(C, 2));
    
    block_matrix<matrix, SymmGroup> op1, op2, op3;
    
    op1 = identity_matrix<matrix>(phys);
    std::cout << "Operator 1:" << std::endl << op1;
    
    {
        matrix tmp(2,2,0.);
        tmp(0,0) = 1.; tmp(1,1) = -1.;
        op2.insert_block(tmp, C, C);
    }
    std::cout << "Operator 2:" << std::endl << op2;
    
    {
        matrix tmp(2,2,0.); tmp(0,1) = 1.;
        op3.insert_block(tmp, C, C);
    }
    std::cout << "Operator 3:" << std::endl << op3;

    all_tests(phys, op1, op2, op3);
}

void test_u1()
{
    typedef U1 SymmGroup;
    maquis::cout << "TESTING U1 SYMMETRY" << std::endl;
    
    Index<SymmGroup> phys;
    phys.insert(std::make_pair(0, 1));
    phys.insert(std::make_pair(1, 2));
    
    block_matrix<matrix, SymmGroup> op1, op2, op3;
    
    op1 = identity_matrix<matrix>(phys);
    std::cout << "Operator 1:" << std::endl << op1;
    
    {
        matrix tmp(1,2,0.); tmp(0,0) = 1.;
        op2.insert_block(tmp, 0, 1);
    }
    std::cout << "Operator 2:" << std::endl << op2;
    
    {
        matrix tmp(2,1,0.); tmp(1,0) = 1.;
        op3.insert_block(tmp, 1, 0);
    }
    std::cout << "Operator 3:" << std::endl << op3;
    
    all_tests(phys, op1, op2, op3);
}


int main() {
    try {
        test_none();
        test_u1();
    } catch (std::exception & e) {
        std::cerr << "Exception:" << std::endl << e.what() << std::endl;
        return 1;
    }
}
