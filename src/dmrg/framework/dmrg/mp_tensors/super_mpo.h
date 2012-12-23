/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2012 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef SUPER_MPO_H
#define SUPER_MPO_H

#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"

#include <boost/function.hpp>

/*
 * Building Super MPS from an MPO object
 *
 * When the code is used for density matrix `rho` evolution, measurements are
 * computed as overlap with a Super MPS.
 * The Super MPS is equivalent to an MPO where the two physical indexes are
 * fused together.
 *
 * Since the MPO doesn't use symmetries for the auxiliary legs, they are mapped
 * to a single block with charge SymmGroup::IdentityCharge.
 *
 * Operators in the MPO are supposed to be in the form:
 *       O_{s1,s2}
 * where s1 is the input state, and s2 the output.
 * (transpose of conventional matrix form)
 * The indexes are fused according to
         s = s1 \otimes adjoin(s2),
 * so that s2 (output) is the most frequent running index.
 *
 * Note: the above mapping is in disagreement with block_to_mpo, but the two
 * functions are anyway not used together.
 */
 
template <class Matrix, class SymmGroup>
MPS<Matrix, SymmGroup> mpo_to_smps(MPO<Matrix, SymmGroup> const& mpo, Index<SymmGroup> const& phys_i)
{
    typedef typename SymmGroup::charge charge;
    typedef boost::unordered_map<size_t,std::pair<charge,size_t> > bond_charge_map;
    
    MPS<Matrix, SymmGroup> mps(mpo.size());
    
    boost::function<charge (charge, charge)> phys_fuse = boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                                                             boost::lambda::_1, -boost::lambda::_2);
    
    Index<SymmGroup> phys2_i = phys_i*adjoin(phys_i);
    ProductBasis<SymmGroup> phys_prod(phys_i, phys_i, phys_fuse);
    Index<SymmGroup> left_i, right_i;
    left_i.insert( std::make_pair(SymmGroup::IdentityCharge, 1) );
    
    bond_charge_map left_map, right_map;
    left_map[0] = std::make_pair(SymmGroup::IdentityCharge, 0);
    
    for (int i=0; i<mpo.size(); ++i) {
        
        ProductBasis<SymmGroup> left_out(phys2_i, left_i);
        boost::unordered_map<charge, size_t> right_sizes;
        
        block_matrix<Matrix, SymmGroup> out_block;
        
        /// run=0 computes the sizes of right blocks
        for (int run=0; run<=1; ++run)
            for (size_t b1=0; b1<mpo[i].row_dim(); ++b1)
            {
                for (size_t b2=0; b2<mpo[i].col_dim(); ++b2)
                {
                    if (!mpo[i].has(b1, b2))
                        continue;
                    
                    /// note: this has to be here, because we don't know if b1 exists
                    charge l_charge; size_t ll;
                    boost::tie(l_charge, ll) = left_map[b1];
                    size_t l_size = left_i[left_i.position(l_charge)].second;
                    
                    block_matrix<Matrix, SymmGroup> const& in_block = mpo[i](b1, b2);
                    for (size_t n=0; n<in_block.n_blocks(); ++n)
                    {
                        charge s1_charge; size_t size1;
                        boost::tie(s1_charge, size1) = in_block.left_basis()[n];
                        charge s2_charge; size_t size2;
                        boost::tie(s2_charge, size2) = in_block.right_basis()[n];
                        
                        charge s_charge = phys_fuse(s1_charge, s2_charge);
                        charge out_l_charge = SymmGroup::fuse(s_charge, l_charge);
                        charge out_r_charge = out_l_charge;
                        
                        if (false && run==0) {
                            std::cout << "s1: " << s1_charge << std::endl;
                            std::cout << "s2: " << s2_charge << std::endl;
                            std::cout << "s:  " << s_charge << std::endl;
                            std::cout << "b1: " << b1 << std::endl;
                            std::cout << "b2: " << b2 << std::endl;
                            std::cout << "l:  " << l_charge << std::endl;
                            std::cout << "r:  " << out_r_charge << std::endl;
                            
                            std::cout << " ------ " << std::endl;
                        }
                        
                        if ( run == 0) {
                            typename bond_charge_map::const_iterator match = right_map.find(b2);
                            if (match == right_map.end()) {
                                right_map[b2] = std::make_pair(out_r_charge, right_sizes[out_r_charge]++);
                            } else
                                assert(match->second.first == out_r_charge);
                            
                            continue;
                        }
                        
                        
                        if (!out_block.has_block(out_l_charge, out_r_charge))
                            out_block.insert_block(Matrix(left_out.size(s_charge,l_charge), right_sizes[out_r_charge], 0.),
                                                   out_l_charge, out_r_charge);
                        
                        size_t phys_offset = phys_prod(s1_charge, s2_charge);
                        size_t left_offset = left_out(s_charge, l_charge);
                        
                        size_t rr = right_map[b2].second;
                        
                        Matrix & out_m = out_block(out_l_charge, out_r_charge);
                        Matrix const& in_m = in_block[n];

                        for (size_t ss2=0; ss2<size2; ++ss2)
                            for (size_t ss1=0; ss1<size1; ++ss1)
                            {
                                size_t ss = ss2 + ss1*size2 + phys_offset;
                                out_m(left_offset + ss*l_size + ll, rr) = in_m(ss1, ss2);
                            }
                    }
                }
            }
        
        right_i = out_block.right_basis();
        
        mps[i] = MPSTensor<Matrix, SymmGroup>(phys2_i, left_i, right_i,
                                              out_block, LeftPaired);
        std::swap(left_i, right_i);
        std::swap(left_map, right_map);
        right_map.clear();
    }
    
    return mps;
}

#endif
