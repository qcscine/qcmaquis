/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2012 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef JOINS_H
#define JOINS_H

#include "dmrg/mp_tensors/mpstensor.h"

template <class Matrix, class SymmGroup>
MPSTensor<Matrix, SymmGroup> join(MPSTensor<Matrix, SymmGroup> const & m1, MPSTensor<Matrix, SymmGroup> const & m2,
                                  boundary_flag_t boundary_f=no_boundary_f)
{
    assert(m1.site_dim() == m2.site_dim());

    Index<SymmGroup> const& phys_i = m1.site_dim();
    
    m1.make_left_paired();
    m2.make_left_paired();
    
    
    MPSTensor<Matrix, SymmGroup> ret;
    // phys_i is untouched
    ret.phys_i = phys_i;
    
    // computing new left_i
    ret.left_i = m1.left_i;
    if (boundary_f != l_boundary_f) {
        for (typename Index<SymmGroup>::const_iterator it = m2.left_i.begin();
             it != m2.left_i.end(); ++it) {
            if (ret.left_i.has(it->first))
                ret.left_i[ret.left_i.position(it->first)].second += it->second;
            else
                ret.left_i.insert(*it);
        }
    }
    
    // computing new right_i
    ret.right_i = m1.right_i;
    if (boundary_f != r_boundary_f) {
        for (typename Index<SymmGroup>::const_iterator it = m2.right_i.begin();
             it != m2.right_i.end(); ++it) {
            if (ret.right_i.has(it->first))
                ret.right_i[ret.right_i.position(it->first)].second += it->second;
            else
                ret.right_i.insert(*it);
        }
    }
    
    ProductBasis<SymmGroup> out_left_pb(phys_i, ret.left_i);
    Index<SymmGroup> const& out_right = ret.right_i;
    
    using std::size_t;
    
    for (size_t t=0; t<2; ++t) // t=0 --> mps1, t=1 --> mps2
    {
        MPSTensor<Matrix, SymmGroup> const & m = (t==0) ? m1 : m2;
        ProductBasis<SymmGroup> in_left(phys_i, m.row_dim());
        
        for (size_t b = 0; b < m.data().n_blocks(); ++b) {
            typename SymmGroup::charge const& sl_charge = m.data().left_basis()[b].first; // phys + left
            typename SymmGroup::charge const& r_charge = m.data().right_basis()[b].first; // right
            size_t out_r_charge_i = out_right.position(r_charge);

            if (!ret.data().has_block(sl_charge, r_charge))
                ret.data().insert_block(Matrix(out_left_pb.size(sl_charge), out_right[out_r_charge_i].second),
                                       sl_charge, r_charge);

            Matrix & nb = ret.data()(sl_charge, r_charge);
            
            size_t in_r_size = m.data().right_basis()[b].second;
            size_t out_r_offset = 0;
            if (t == 1 && boundary_f != r_boundary_f)
                out_r_offset += m1.col_dim().size_of_block(r_charge, true);
            
            for (size_t s=0; s<phys_i.size(); ++s) {
                typename SymmGroup::charge const& s_charge = phys_i[s].first;
                typename SymmGroup::charge l_charge = SymmGroup::fuse(sl_charge, -s_charge); // left
                
                if (!m.row_dim().has(l_charge))
                    continue;
                
                size_t in_l_size = m.row_dim().size_of_block(l_charge, true);
                size_t in_l_offset = in_left(s_charge, l_charge);

                size_t out_l_size = ret.row_dim().size_of_block(l_charge, true);
                size_t out_l_offset = out_left_pb(s_charge, l_charge);
                if (t == 1 && boundary_f != l_boundary_f)
                    out_l_offset += m1.row_dim().size_of_block(l_charge, true);
                
                for (size_t ss=0; ss<phys_i[s].second; ++ss) {
#ifdef USE_AMBIENT
    assert(false);
#else
                    maquis::dmrg::detail::copy2d(nb, out_l_offset+ss*out_l_size, out_r_offset, m.data()[b], in_l_offset+ss*in_l_size, 0,
                                                 in_l_size, in_r_size);
#endif
                }
            }
        }
    }
    
    // check right_pairing
    assert(ret.right_i == ret.data_.right_basis());
    
    return ret;
}


#endif
