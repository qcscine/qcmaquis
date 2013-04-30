/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

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
            charge new_left = phys_fuse(A.left_basis()[i].first, B.left_basis()[j].first);
            charge new_right = phys_fuse(A.right_basis()[i].first, B.right_basis()[j].first);
            
            Matrix tmp(pb_left.size(new_left), pb_right.size(new_right), 0);
            
            maquis::dmrg::detail::op_kron(tmp, B[j], A[i],
                                          pb_left(A.left_basis()[i].first, B.left_basis()[j].first),
                                          pb_right(A.right_basis()[i].first, B.right_basis()[j].first),
                                          A.left_basis()[i].second, B.left_basis()[j].second,
                                          A.right_basis()[i].second, B.right_basis()[j].second);
            
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
            charge new_left  = group(A.left_basis()[i].first, -B.left_basis()[j].first);
            charge new_right = group(A.right_basis()[i].first, -B.right_basis()[j].first);
            
            Matrix tmp(phys_rho.size_of_block(new_left), phys_rho.size_of_block(new_right), 0);
            
            maquis::dmrg::detail::op_kron(tmp, B[j], A[i], 0, 0,
                                          A.left_basis()[i].second, B.left_basis()[j].second,
                                          A.right_basis()[i].second, B.right_basis()[j].second);
            
            C.match_and_add_block(tmp, new_left, new_right);
        }
    }
}

#endif
