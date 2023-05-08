/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef MPS_MPO_DETAIL_H
#define MPS_MPO_DETAIL_H

namespace mps_mpo_detail {

template<class Matrix, class SymmGroup>
Boundary<Matrix, SymmGroup>
mixed_left_boundary(MPS<Matrix, SymmGroup> const & bra, MPS<Matrix, SymmGroup> const & ket)
{
    assert(ket.length() == bra.length());
    Index<SymmGroup> i = ket[0].row_dim();
    Index<SymmGroup> j = bra[0].row_dim();
    Boundary<Matrix, SymmGroup> ret(i, j, 1);
    for(typename Index<SymmGroup>::basis_iterator it1 = i.basis_begin(); !it1.end(); ++it1)
        for(typename Index<SymmGroup>::basis_iterator it2 = j.basis_begin(); !it2.end(); ++it2)
            ret[0](*it1, *it2) = 1;
    return ret;
}

template<class Matrix, class SymmGroup>
Boundary<Matrix, SymmGroup>
mixed_right_boundary(MPS<Matrix, SymmGroup> const & bra, MPS<Matrix, SymmGroup> const & ket)
{
    assert(ket.length() == bra.length());
    Index<SymmGroup> j = ket[ket.length()-1].col_dim();
    Index<SymmGroup> i = bra[bra.length()-1].col_dim();
    Boundary<Matrix, SymmGroup> ret(i, j, 1);
    for(typename Index<SymmGroup>::basis_iterator it1 = i.basis_begin(); !it1.end(); ++it1)
        for(typename Index<SymmGroup>::basis_iterator it2 = j.basis_begin(); !it2.end(); ++it2)
            ret[0](*it1, *it2) = 1;
    return ret;
}

} // mps_mpo_detail

#endif
