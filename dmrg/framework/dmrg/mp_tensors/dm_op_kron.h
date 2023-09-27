/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef MAQUIS_DMRG_DM_OP_KRON_H
#define MAQUIS_DMRG_DM_OP_KRON_H

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/grouped_symmetry.h"

template<class Matrix, class SymmGroup>
void dm_kron(Index<SymmGroup> const & phys,
             block_matrix<Matrix, SymmGroup> const & A,
             block_matrix<Matrix, SymmGroup> const & B,
             block_matrix<Matrix, SymmGroup> & C)
{
    C = block_matrix<Matrix, SymmGroup>();
    
    Index<SymmGroup> const & left_basis = phys;
    Index<SymmGroup> const & right_basis = phys;
    
    typedef typename SymmGroup::charge charge;
    boost::function<charge (charge, charge)> phys_fuse = boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                                                             boost::lambda::_1, -boost::lambda::_2);
    ProductBasis<SymmGroup> pb_left(left_basis, left_basis, phys_fuse);
    ProductBasis<SymmGroup> const& pb_right = pb_left;
    
    for (int i=0; i<A.n_blocks(); ++i) {
        for (int j=0; j<B.n_blocks(); ++j) {
            charge new_left = phys_fuse(A.basis().left_charge(i), B.basis().left_charge(j));
            charge new_right = phys_fuse(A.basis().right_charge(i), B.basis().right_charge(j));
            
            Matrix tmp(pb_left.size(new_left), pb_right.size(new_right), 0);
            
            maquis::dmrg::detail::op_kron(tmp, B[j], A[i],
                                          pb_left(A.basis().left_charge(i), B.basis().left_charge(j)),
                                          pb_right(A.basis().right_charge(i), B.basis().right_charge(j)),
                                          A.basis().left_size(i), B.basis().left_size(j),
                                          A.basis().right_size(i), B.basis().right_size(j));
            
            C.match_and_add_block(tmp, new_left, new_right);
        }
    }
}

template<class Matrix, class SymmGroup>
void dm_group_kron(Index<SymmGroup> const & phys_psi,
                   block_matrix<Matrix, SymmGroup> const & A,
                   block_matrix<Matrix, SymmGroup> const & B,
                   block_matrix<Matrix, typename grouped_symmetry<SymmGroup>::type> & C)
{
    typedef typename grouped_symmetry<SymmGroup>::type OutSymm;
    C = block_matrix<Matrix, OutSymm>();
    
    Index<SymmGroup> const & left_basis = phys_psi;
    Index<SymmGroup> const & right_basis = phys_psi;
    
    typedef typename OutSymm::charge charge;
    Index<OutSymm> phys_rho = group(phys_psi, adjoin(phys_psi));
    
    for (int i=0; i<A.n_blocks(); ++i) {
        for (int j=0; j<B.n_blocks(); ++j) {
            charge new_left  = group(A.basis().left_charge(i), -B.basis().left_charge(j));
            charge new_right = group(A.basis().right_charge(i), -B.basis().right_charge(j));
            
            Matrix tmp(phys_rho.size_of_block(new_left), phys_rho.size_of_block(new_right), 0);
            
            maquis::dmrg::detail::op_kron(tmp, B[j], A[i], 0, 0,
                                          A.basis().left_size(i), B.basis().left_size(j),
                                          A.basis().right_size(i), B.basis().right_size(j));
            
            C.match_and_add_block(tmp, new_left, new_right);
        }
    }
}

#endif
