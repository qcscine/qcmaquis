/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2021 Institute for Theoretical Physics, ETH Zurich
 *               2021 by Alberto Baiardi <abaiardi@ethz.ch>
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