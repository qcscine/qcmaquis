/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef MPS_REFLECT_H
#define MPS_REFLECT_H

#include "dmrg/mp_tensors/mps.h"

template <class Matrix, class SymmGroup>
MPS<Matrix, SymmGroup> reflect(MPS<Matrix, SymmGroup> const& mps_in)
{
    typedef std::size_t size_type;
    typedef typename SymmGroup::charge charge;
    const size_type L = mps_in.length();
    const charge final_charge = mps_in.col_dim(L-1)[0].first;
    
    MPS<Matrix, SymmGroup> ret(L);
    for (size_type p=0; p<L; ++p) {
        
        mps_in[p].make_left_paired();
        
        
        Index<SymmGroup> const& phys_i  = mps_in[p].site_dim();
        Index<SymmGroup> const& left_i  = mps_in[p].row_dim();
        Index<SymmGroup> const& right_i = mps_in[p].col_dim();
        
        ProductBasis<SymmGroup> in_left(phys_i, left_i);
        
        Index<SymmGroup> out_left_i(right_i), out_right_i(left_i);
        for (size_type i=0; i<out_left_i.size(); ++i)
            out_left_i[i].first = SymmGroup::fuse(final_charge, -out_left_i[i].first);
        out_left_i.sort();
        for (size_type i=0; i<out_right_i.size(); ++i)
            out_right_i[i].first = SymmGroup::fuse(final_charge, -out_right_i[i].first);
        out_right_i.sort();
        
        ProductBasis<SymmGroup> out_left(phys_i, out_left_i);
        
        block_matrix<Matrix, SymmGroup> const& m = mps_in[p].data();
        block_matrix<Matrix, SymmGroup> mnew;
        
        for (size_t block = 0; block < m.n_blocks(); ++block) {
            const size_type r = right_i.position(m.basis().right_charge(block));
            if(r == right_i.size()) continue;
            const charge in_l_charge = m.basis().left_charge(block);
            const charge in_r_charge = right_i[r].first;
            
            const size_type out_l = out_left_i.position(SymmGroup::fuse(final_charge, -right_i[r].first));
            
            for (size_t s = 0; s < phys_i.size(); ++s) {
                const size_type l = left_i.position(SymmGroup::fuse(in_l_charge, -phys_i[s].first));
                if(l == left_i.size()) continue;
                
                const size_type out_r = out_right_i.position(SymmGroup::fuse(final_charge, -left_i[l].first));
                
                const charge out_l_charge = SymmGroup::fuse(phys_i[s].first, out_left_i[out_l].first);
                const charge out_r_charge = out_right_i[out_r].first;
                
                if(!mnew.has_block(out_l_charge, out_r_charge))
                    mnew.insert_block(Matrix(out_left.size(out_l_charge), out_right_i[out_r].second, 0),
                                      out_l_charge, out_r_charge);
                
                const size_type in_left_offset = in_left(phys_i[s].first, left_i[l].first);
                const size_type out_left_offset = out_left(phys_i[s].first, out_left_i[out_l].first);
                
                Matrix const& in_block = m[block];
                Matrix & out_block = mnew(out_l_charge, out_r_charge);
                
                for (size_type ss=0; ss<phys_i[s].second; ++ss) {
                    for (size_type ll=0; ll<left_i[l].second; ++ll) {
                        for (size_type rr=0; rr<right_i[r].second; ++rr) {
                            
                            const size_type ii_in  = ss*left_i[l].second + in_left_offset + ll;
                            const size_type jj_in  = rr;
                            const size_type ii_out = ss*out_left_i[out_l].second + out_left_offset + rr;
                            const size_type jj_out = ll;
                            
                            out_block(ii_out, jj_out) = in_block(ii_in, jj_in);
                        }
                    }
                }
            }
        }
        
        ret[L-1 - p] = MPSTensor<Matrix, SymmGroup>(phys_i, out_left_i, out_right_i, mnew, LeftPaired);
        
    }
    
    return ret;
    
    
}

#endif
