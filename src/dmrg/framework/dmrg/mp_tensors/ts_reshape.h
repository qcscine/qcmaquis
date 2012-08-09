/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#include "twositetensor.h"
#include "dmrg/mp_tensors/mpstensor.h"

#include "dmrg/mp_tensors/reshapes.h"
#include "dmrg/mp_tensors/generic_reshape.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"

#include "dmrg/utils/random.hpp"
#include <alps/numeric/real.hpp>

#include <vector>
#include <utility>

namespace ts_reshape {
    
    template<class Matrix, class SymmGroup>
    void reshape_both_to_left(Index<SymmGroup> physical_i,
                              Index<SymmGroup> left_i,
                              Index<SymmGroup> right_i,
                              block_matrix<Matrix, SymmGroup> const & m1,
                              block_matrix<Matrix, SymmGroup> & m2)
    {   
        
        m2 = block_matrix<Matrix, SymmGroup>();
        
        typedef std::size_t size_t;
        typedef typename SymmGroup::charge charge;
        
        Index<SymmGroup> phys2_i = physical_i*physical_i;
        ProductBasis<SymmGroup> phys_pb(physical_i, physical_i);
        ProductBasis<SymmGroup> in_left(physical_i, left_i);
        ProductBasis<SymmGroup> in_right(physical_i, right_i,
                                         boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                                             -boost::lambda::_1, boost::lambda::_2));
        ProductBasis<SymmGroup> out_left(phys2_i, left_i);
        
       
        for (size_t block = 0; block < m1.n_blocks(); ++block)
        {
            for (size_t s1 = 0; s1 < physical_i.size(); ++s1)
                for (size_t s2 = 0; s2 < physical_i.size(); ++s2)
                {
                    size_t l = left_i.position(SymmGroup::fuse(m1.left_basis()[block].first,
                                                              -physical_i[s1].first));
                    if(l == left_i.size()) continue;
                    size_t r = right_i.position(SymmGroup::fuse(m1.right_basis()[block].first,
                                                                physical_i[s2].first));
                    if(r == right_i.size()) continue;
                    
                    {                            
                        charge s_charge = SymmGroup::fuse(physical_i[s1].first, physical_i[s2].first);
                        size_t s_out = phys2_i.position(s_charge);
                        
                        charge in_l_charge = SymmGroup::fuse(physical_i[s1].first, left_i[l].first);
                        charge in_r_charge = SymmGroup::fuse(-physical_i[s2].first, right_i[r].first);
                        charge out_l_charge = SymmGroup::fuse(s_charge, left_i[l].first);
                        charge out_r_charge = right_i[r].first;
                        
                        if (! m1.has_block(in_l_charge, in_r_charge) ) continue;
                        
                        // Why is this supposed to work?
                        if (! m2.has_block(out_l_charge, out_r_charge) )
                            m2.insert_block(new Matrix(out_left.size(s_charge, left_i[l].first), right_i[r].second, 0),
                                            out_l_charge, out_r_charge);
                        
                        maquis::dmrg::detail::reshape_b2l( m2(out_l_charge, out_r_charge), m1(in_l_charge, in_r_charge), 
                                                           in_left(physical_i[s1].first, left_i[l].first), in_right(physical_i[s2].first, right_i[r].first),
                                                           out_left(s_charge, left_i[l].first), phys_pb(physical_i[s1].first, physical_i[s2].first),
                                                           physical_i[s1].second, physical_i[s2].second, left_i[l].second, right_i[r].second );
                    }
                }
        }
        
    }
    
    /* 
    template<class Matrix, class SymmGroup>
    void reshape_right_to_left(Index<SymmGroup> physical_i,
                               Index<SymmGroup> left_i,
                               Index<SymmGroup> right_i,
                               block_matrix<Matrix, SymmGroup> const & m1,
                               block_matrix<Matrix, SymmGroup> & m2)
    {   }
    */
    
    template<class Matrix, class SymmGroup>
    void reshape_left_to_both(Index<SymmGroup> physical_i,
                              Index<SymmGroup> left_i,
                              Index<SymmGroup> right_i,
                              block_matrix<Matrix, SymmGroup> const & m1,
                              block_matrix<Matrix, SymmGroup> & m2)
    {   
        
        m2 = block_matrix<Matrix, SymmGroup>();
        
        typedef std::size_t size_t;
        typedef typename SymmGroup::charge charge;
        
        Index<SymmGroup> phys2_i = physical_i*physical_i;
        ProductBasis<SymmGroup> phys_pb(physical_i, physical_i);
        ProductBasis<SymmGroup> in_left(phys2_i, left_i);
        
        ProductBasis<SymmGroup> out_right(physical_i, right_i,
                                          boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                                              -boost::lambda::_1, boost::lambda::_2));
        ProductBasis<SymmGroup> out_left(physical_i, left_i);
       
        for (size_t block = 0; block < m1.n_blocks(); ++block)
        {
            for (size_t s1 = 0; s1 < physical_i.size(); ++s1)
                for (size_t s2 = 0; s2 < physical_i.size(); ++s2)
                {
                    charge s_charge = SymmGroup::fuse(physical_i[s1].first, physical_i[s2].first);
                    size_t s_out = phys2_i.position(s_charge);
                    
                    size_t r = right_i.position(m1.right_basis()[block].first);
                    if(r == right_i.size()) continue;
                    size_t l = left_i.position(SymmGroup::fuse(m1.left_basis()[block].first, -s_charge));
                    if(l == left_i.size()) continue;

                    {
                        charge out_l_charge = SymmGroup::fuse(physical_i[s1].first, left_i[l].first);
                        charge out_r_charge = SymmGroup::fuse(-physical_i[s2].first, right_i[r].first);
                        charge in_l_charge = SymmGroup::fuse(s_charge, left_i[l].first);
                        charge in_r_charge = right_i[r].first;
                        
                        if (! m1.has_block(in_l_charge, in_r_charge) ) continue;
                        
                        if (! m2.has_block(out_l_charge, out_r_charge) )
                            m2.insert_block(new Matrix(out_left.size(physical_i[s1].first, left_i[l].first),
                                                   out_right.size(-physical_i[s2].first, right_i[r].first), 0),
                                            out_l_charge, out_r_charge);
                        
                        maquis::dmrg::detail::reshape_l2b( m2(out_l_charge, out_r_charge), m1(in_l_charge, in_r_charge),
                                                           in_left(s_charge, left_i[l].first), phys_pb(physical_i[s1].first, physical_i[s2].first),
                                                           out_left(physical_i[s1].first, left_i[l].first), out_right(physical_i[s2].first, right_i[r].first),
                                                           physical_i[s1].second, physical_i[s2].second, left_i[l].second, right_i[r].second );
                    }
                }
        }
        
    }
    
    template<class Matrix, class SymmGroup>
    void reshape_right_to_both(Index<SymmGroup> physical_i,
                               Index<SymmGroup> left_i,
                               Index<SymmGroup> right_i,
                               block_matrix<Matrix, SymmGroup> const & m1,
                               block_matrix<Matrix, SymmGroup> & m2)
    {
        
        m2 = block_matrix<Matrix, SymmGroup>();
        
        typedef std::size_t size_t;
        typedef typename SymmGroup::charge charge;
        
        Index<SymmGroup> phys2_i = physical_i*physical_i;
        ProductBasis<SymmGroup> phys_pb(physical_i, physical_i);
        ProductBasis<SymmGroup> in_right(phys2_i, right_i,
                                         boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                                             -boost::lambda::_1, boost::lambda::_2));
        
        ProductBasis<SymmGroup> out_right(physical_i, right_i,
                                          boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                                              -boost::lambda::_1, boost::lambda::_2));
        ProductBasis<SymmGroup> out_left(physical_i, left_i);
        
        for (size_t block = 0; block < m1.n_blocks(); ++block)
        {
            for (size_t s1 = 0; s1 < physical_i.size(); ++s1)
                for (size_t s2 = 0; s2 < physical_i.size(); ++s2)
                {
                    charge s_charge = SymmGroup::fuse(physical_i[s1].first, physical_i[s2].first);
                    size_t s_out = phys2_i.position(s_charge);
                  
                    size_t l = left_i.position(m1.left_basis()[block].first);
                    if(l == left_i.size()) continue;
                    size_t r = right_i.position(SymmGroup::fuse(m1.right_basis()[block].first, s_charge));
                    if(r == right_i.size()) continue;

                    {
                        charge out_l_charge = SymmGroup::fuse(physical_i[s1].first, left_i[l].first);
                        charge out_r_charge = SymmGroup::fuse(-physical_i[s2].first, right_i[r].first);
                        charge in_l_charge = left_i[l].first;
                        charge in_r_charge = SymmGroup::fuse(-s_charge, right_i[r].first);
                        
                        if (! m1.has_block(in_l_charge, in_r_charge) ) continue;
                        
                        if (! m2.has_block(out_l_charge, out_r_charge) )
                            m2.insert_block(new Matrix(out_left.size(physical_i[s1].first, left_i[l].first),
                                                   out_right.size(-physical_i[s2].first, right_i[r].first), 0),
                                            out_l_charge, out_r_charge);
                        
                        size_t in_right_offset  = in_right  (s_charge            , right_i[r].first);
                        size_t out_right_offset = out_right (physical_i[s2].first, right_i[r].first);
                        size_t out_left_offset  = out_left  (physical_i[s1].first, left_i[l].first);
                        size_t in_phys_offset   = phys_pb   (physical_i[s1].first, physical_i[s2].first);
                        
                        Matrix const & in_block = m1(in_l_charge, in_r_charge);
                        Matrix & out_block = m2(out_l_charge, out_r_charge);
                        
                        for (size_t ss1 = 0; ss1 < physical_i[s1].second; ++ss1)
                            for (size_t ss2 = 0; ss2 < physical_i[s2].second; ++ss2)
                            {
                                size_t ss_out = in_phys_offset + ss1*physical_i[s2].second + ss2;
                                for (size_t rr = 0; rr < right_i[r].second; ++rr)
                                    for (size_t ll = 0; ll < left_i[l].second; ++ll)
                                        out_block(out_left_offset + ss1*left_i[l].second + ll, out_right_offset + ss2*right_i[r].second + rr) = in_block(ll, in_right_offset + ss_out*right_i[r].second + rr);
                            }
                        
                    }
                }
        }
    }
    
    /*
    template<class Matrix, class SymmGroup>
    void reshape_left_to_right(Index<SymmGroup> physical_i,
                               Index<SymmGroup> left_i,
                               Index<SymmGroup> right_i,
                               block_matrix<Matrix, SymmGroup> const & m1,
                               block_matrix<Matrix, SymmGroup> & m2)
    {   }
    */
    
    template<class Matrix, class SymmGroup>
    void reshape_both_to_right(Index<SymmGroup> physical_i,
                               Index<SymmGroup> left_i,
                               Index<SymmGroup> right_i,
                               block_matrix<Matrix, SymmGroup> const & m1,
                               block_matrix<Matrix, SymmGroup> & m2)
    {   
        
        m2 = block_matrix<Matrix, SymmGroup>();
        
        typedef std::size_t size_t;
        typedef typename SymmGroup::charge charge;
        
        Index<SymmGroup> phys2_i = physical_i*physical_i;
        ProductBasis<SymmGroup> phys_pb(physical_i, physical_i);
        ProductBasis<SymmGroup> in_left(physical_i, left_i);
        ProductBasis<SymmGroup> in_right(physical_i, right_i,
                                         boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                                             -boost::lambda::_1, boost::lambda::_2));
        ProductBasis<SymmGroup> out_right(phys2_i, right_i,
                                          boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                                              -boost::lambda::_1, boost::lambda::_2));
       
        for (size_t block = 0; block < m1.n_blocks(); ++block)
        {
            for (size_t s1 = 0; s1 < physical_i.size(); ++s1)
                for (size_t s2 = 0; s2 < physical_i.size(); ++s2)
                {
                    size_t l = left_i.position(SymmGroup::fuse(m1.left_basis()[block].first,
                                                               -physical_i[s1].first));
                    if(l == left_i.size()) continue;
                    size_t r = right_i.position(SymmGroup::fuse(m1.right_basis()[block].first,
                                                                physical_i[s2].first));
                    if(r == right_i.size()) continue;

                    {
                        charge s_charge = SymmGroup::fuse(physical_i[s1].first, physical_i[s2].first);
                        size_t s_out = phys2_i.position(s_charge);
                        
                        charge in_l_charge = SymmGroup::fuse(physical_i[s1].first, left_i[l].first);
                        charge in_r_charge = SymmGroup::fuse(-physical_i[s2].first, right_i[r].first);
                        charge out_l_charge = left_i[l].first;
                        charge out_r_charge = SymmGroup::fuse(-s_charge, right_i[r].first);
                        
                        if (! m1.has_block(in_l_charge, in_r_charge) ) continue;
                        
                        if (! m2.has_block(out_l_charge, out_r_charge) )
                            m2.insert_block(new Matrix(left_i[l].second, out_right.size(-s_charge, right_i[r].first), 0),
                                            out_l_charge, out_r_charge);
                        
                        size_t in_left_offset = in_left(physical_i[s1].first, left_i[l].first);
                        size_t in_right_offset = in_right(physical_i[s2].first, right_i[r].first);
                        size_t out_right_offset = out_right(s_charge, right_i[r].first);
                        size_t out_phys_offset = phys_pb(physical_i[s1].first, physical_i[s2].first);
                        
                        Matrix const & in_block = m1(in_l_charge, in_r_charge);
                        Matrix & out_block = m2(out_l_charge, out_r_charge);
                        
                        for (size_t ss1 = 0; ss1 < physical_i[s1].second; ++ss1)
                            for (size_t ss2 = 0; ss2 < physical_i[s2].second; ++ss2)
                            {
                                size_t ss_out = out_phys_offset + ss1*physical_i[s2].second + ss2;
                                for (size_t rr = 0; rr < right_i[r].second; ++rr)
                                    for (size_t ll = 0; ll < left_i[l].second; ++ll)
                                        out_block(ll, out_right_offset + ss_out*right_i[r].second + rr) = in_block(in_left_offset + ss1*left_i[l].second+ll, in_right_offset + ss2*right_i[r].second+rr);
                            }
                    }
                }
        }
        
    }
    
} // namespace ts_reshape
