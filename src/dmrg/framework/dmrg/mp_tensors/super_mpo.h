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
#include "dmrg/block_matrix/grouped_symmetry.h"

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
    typedef typename MPOTensor<Matrix, SymmGroup>::row_proxy row_proxy;
    
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
                //for (size_t b2=0; b2<mpo[i].col_dim(); ++b2)
                row_proxy row_b1 = mpo[i].row(b1);
                for (typename row_proxy::const_iterator it = row_b1.begin(); it != row_b1.end(); ++it)
                {
                    size_t b2 = it.index(); 
                    
                    /// note: this has to be here, because we don't know if b1 exists
                    charge l_charge; size_t ll;
                    boost::tie(l_charge, ll) = left_map[b1];
                    size_t l_size = left_i[left_i.position(l_charge)].second;
                    
                    typename Matrix::value_type scale = mpo[i].at(b1, b2).scale;
                    block_matrix<Matrix, SymmGroup> const& in_block = mpo[i].at(b1, b2).op;
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
                            maquis::cout << "s1: " << s1_charge << std::endl;
                            maquis::cout << "s2: " << s2_charge << std::endl;
                            maquis::cout << "s:  " << s_charge << std::endl;
                            maquis::cout << "b1: " << b1 << std::endl;
                            maquis::cout << "b2: " << b2 << std::endl;
                            maquis::cout << "l:  " << l_charge << std::endl;
                            maquis::cout << "r:  " << out_r_charge << std::endl;
                            
                            maquis::cout << " ------ " << std::endl;
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
                                // TODO: Sebastian thinks, this is correct but should be checked
                                out_m(left_offset + ss*l_size + ll, rr) = in_m(ss1, ss2) * scale;
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
 * The indexes are fused/grouped according to
 s = <s1, adjoin(s2)>,
 * so that s2 (output) is the most frequent running index.
 */

template <class Matrix, class InSymm>
MPS<Matrix, typename grouped_symmetry<InSymm>::type> mpo_to_smps_group(MPO<Matrix, InSymm> const& mpo, Index<InSymm> const& phys_i,
                                                                       std::vector<Index<typename grouped_symmetry<InSymm>::type> > const& allowed)
{
    typedef typename grouped_symmetry<InSymm>::type OutSymm;
    typedef typename InSymm::charge in_charge;
    typedef typename OutSymm::charge out_charge;
    typedef boost::unordered_map<size_t,std::pair<out_charge,size_t> > bond_charge_map;
    typedef typename MPOTensor<Matrix, InSymm>::row_proxy row_proxy;
    
    MPS<Matrix, OutSymm> mps(mpo.size());
    
    boost::function<out_charge (in_charge, in_charge)> phys_group = boost::lambda::bind(static_cast<out_charge(*)(in_charge, in_charge)>(group),
                                                                                        boost::lambda::_1, -boost::lambda::_2);
    
    Index<OutSymm> phys2_i = group(phys_i, adjoin(phys_i));
    Index<OutSymm> left_i, right_i;
    left_i.insert( std::make_pair(OutSymm::IdentityCharge, 1) );

    for (int i=0; i<mpo.size(); ++i) {
        ProductBasis<OutSymm> left_out(phys2_i, left_i);
        
        block_matrix<Matrix, OutSymm> out_block;
            //for (size_t b1=0; b1<mpo[i].row_dim(); ++b1)
            //{
            //    for (size_t b2=0; b2<mpo[i].col_dim(); ++b2)
            //    {
            //        if (!mpo[i].has(b1, b2))
            //            continue;
            
            for (size_t b1=0; b1<mpo[i].row_dim(); ++b1)
            {
                row_proxy row_b1 = mpo[i].row(b1);
                for (typename row_proxy::const_iterator it = row_b1.begin(); it != row_b1.end(); ++it)
                {
                    size_t b2 = it.index(); 

                    for (size_t l=0; l<left_i.size(); ++l)
                    {
                        out_charge l_charge = left_i[l].first;
                        size_t     l_size   = left_i[l].second;
                        size_t     ll       = b1;
                        
                        typename Matrix::value_type scale = mpo[i].at(b1, b2).scale;
                        block_matrix<Matrix, InSymm> const& in_block = mpo[i].at(b1, b2).op;
                        for (size_t n=0; n<in_block.n_blocks(); ++n)
                        {
                            in_charge s1_charge; size_t size1;
                            boost::tie(s1_charge, size1) = in_block.left_basis()[n];
                            in_charge s2_charge; size_t size2;
                            boost::tie(s2_charge, size2) = in_block.right_basis()[n];
                            
                            out_charge s_charge = phys_group(s1_charge, s2_charge);
                            out_charge out_l_charge = OutSymm::fuse(s_charge, l_charge);
                            out_charge out_r_charge = out_l_charge;
                            
                            if (! allowed[i+1].has(out_r_charge) )
                                continue;
                            
                            if (!out_block.has_block(out_l_charge, out_r_charge))
                                out_block.insert_block(Matrix(left_out.size(s_charge,l_charge), mpo[i].col_dim(), 0.),
                                                       out_l_charge, out_r_charge);
                            
                            size_t phys_offset = 0;
                            size_t left_offset = left_out(s_charge, l_charge);
                            
                            size_t rr = b2;
                            
                            Matrix & out_m = out_block(out_l_charge, out_r_charge);
                            Matrix const& in_m = in_block[n];
                            
                            for (size_t ss2=0; ss2<size2; ++ss2)
                                for (size_t ss1=0; ss1<size1; ++ss1)
                                {
                                    size_t ss = ss2 + ss1*size2 + phys_offset;
                                    // TODO: Sebastian thinks, this is correct but should be checked
                                    out_m(left_offset + ss*l_size + ll, rr) = in_m(ss1, ss2) * scale;
                                }
                        }
                    }
                }
            }
        
        right_i = out_block.right_basis();
        
        mps[i] = MPSTensor<Matrix, OutSymm>(phys2_i, left_i, right_i,
                                            out_block, LeftPaired);
        std::swap(left_i, right_i);
    }
    
    return mps;
}

#endif
