/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef RESHAPE_H
#define RESHAPE_H

#include <map>

#include "block_matrix/indexing.h"
#include "utils/ambient_assert.h"

template<class Matrix, class SymmGroup>
void reshape_left_to_right(Index<SymmGroup> physical_i,
                           Index<SymmGroup> left_i,
                           Index<SymmGroup> right_i,
                           block_matrix<Matrix, SymmGroup> const & m1,
                           block_matrix<Matrix, SymmGroup> & m2)
{   
    static Timer timer("reshape_left_to_right");
    timer.begin();
    
    using std::size_t;
    
    m2 = block_matrix<Matrix, SymmGroup>();
    
    typedef std::size_t size_t;
    typedef typename SymmGroup::charge charge;
    
    for (int run = 0; run < 2; ++run) {
        ProductBasis<SymmGroup> in_left(physical_i, left_i);
        ProductBasis<SymmGroup> out_right(physical_i, right_i,
                                          boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                                              -boost::lambda::_1, boost::lambda::_2));
        
        if (run == 1)
            m2.allocate_blocks();
       
        for (size_t block = 0; block < m1.n_blocks(); ++block)
        {
            for (size_t s = 0; s < physical_i.size(); ++s)
            {
                size_t l = left_i.position(SymmGroup::fuse(m1.left_basis()[block].first,
                                                           -physical_i[s].first));
                size_t r = right_i.position(m1.right_basis()[block].first);
                
                if (l == left_i.size())
                    continue;
                if (r == right_i.size())
                    continue;
                
                {
                    bool pretend = (run == 0);
                    
                    charge in_l_charge = SymmGroup::fuse(physical_i[s].first, left_i[l].first);
                    charge in_r_charge = right_i[r].first;
                    charge out_l_charge = left_i[l].first;
                    charge out_r_charge = SymmGroup::fuse(-physical_i[s].first, right_i[r].first);

                    if (! m1.has_block(in_l_charge, in_r_charge) )
                        continue;
                    
                    size_t in_left_offset = in_left(physical_i[s].first, left_i[l].first);
                    size_t out_right_offset = out_right(physical_i[s].first, right_i[r].first);
                    
                    if (!pretend){
                        Matrix const & in_block = m1(in_l_charge, in_r_charge);
                        Matrix & out_block = m2(out_l_charge, out_r_charge);
                        
                        /* optimize me */
                        #ifdef MPI_PARALLEL
                        ambient::push(ambient::reshape_l2r_l_kernel, ambient::reshape_l2r_c_kernel, in_block, out_block, 
                                      in_left_offset, out_right_offset, physical_i[s].second, left_i[l].second, right_i[r].second);
                        #else
                        for (size_t ss = 0; ss < physical_i[s].second; ++ss)
                            for (size_t rr = 0; rr < right_i[r].second; ++rr)
                                for (size_t ll = 0; ll < left_i[l].second; ++ll)
                                    out_block(ll, out_right_offset + ss*right_i[r].second+rr) = in_block(in_left_offset + ss*left_i[l].second+ll, rr);
                        #endif
                    }

                    if (pretend)
                        m2.reserve(out_l_charge, out_r_charge,
                                   left_i[l].second, out_right_offset + physical_i[s].second * right_i[r].second);
                }
            }
        }
    }
    
    timer.end();
    
//    assert(m2.left_basis() == left_i);
}

/*
 pseudo-code kernel
 
 std::map<std::pair<charge, charge>, std::list<size_t> > block_requirements;
 
 for (s = 0; s < physical_i.size(); ++s)
	 for (l = 0; ...)
		 for (r = 0; ...) {
             charge in_l_charge = SymmGroup::fuse(physical_i[s].first, left_i[l].first);
             charge in_r_charge = right_i[r].first;
 
			 charge out_l_charge = left_i[l].first;
			 charge out_r_charge = SymmGroup::fuse(-physical_i[s].first, right_i[r].first);
 			
 			 block_requirements[make_pair(out_l_charge, out_r_charge)].push_back(m1.position(in_l_charge, in_r_charge);
 }
 
 A: logistics kernel
 
 // on the node that creates (l,r) in the final block_matrix,
 // get block_requirements[(l,r)] from other nodes
 
 B: compute kernel
 
 */
 

template<class Matrix, class SymmGroup>
void reshape_right_to_left(Index<SymmGroup> physical_i,
                           Index<SymmGroup> left_i,
                           Index<SymmGroup> right_i,
                           block_matrix<Matrix, SymmGroup> const & m1,
                           block_matrix<Matrix, SymmGroup> & m2)
{
//    assert(m1.left_basis() == left_i);
    
    static Timer timer("reshape_right_to_left");
    timer.begin();
    
    using std::size_t;
    
    m2 = block_matrix<Matrix, SymmGroup>();
    
    typedef std::size_t size_t;
    typedef typename SymmGroup::charge charge;
    
    for (int run = 0; run < 2; ++run) {
        ProductBasis<SymmGroup> in_right(physical_i, right_i,
                                         boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                                             -boost::lambda::_1, boost::lambda::_2));
        ProductBasis<SymmGroup> out_left(physical_i, left_i);
        
        if (run == 1)
            m2.allocate_blocks();
    
        for (size_t block = 0; block < m1.n_blocks(); ++block)
        {
            for (size_t s = 0; s < physical_i.size(); ++s)
            {
                size_t l = left_i.position(m1.left_basis()[block].first);
                size_t r = right_i.position(SymmGroup::fuse(m1.right_basis()[block].first,
                                                            physical_i[s].first));
                if(l == left_i.size()) continue;
                if(r == right_i.size()) continue;
                {
                    bool pretend = (run == 0);
                    
                    charge in_l_charge = left_i[l].first;
                    charge in_r_charge = SymmGroup::fuse(-physical_i[s].first, right_i[r].first);
                    charge out_l_charge = SymmGroup::fuse(physical_i[s].first, left_i[l].first);
                    charge out_r_charge = right_i[r].first;
                    
                    if (! m1.has_block(in_l_charge, in_r_charge) )
                        continue;
                    
                    size_t in_right_offset = in_right(physical_i[s].first, right_i[r].first);
                    size_t out_left_offset = out_left(physical_i[s].first, left_i[l].first);
                    
                    if (!pretend) {
                        Matrix const & in_block = m1(in_l_charge, in_r_charge);
                        Matrix & out_block = m2(out_l_charge, out_r_charge);
                        
                        /* optimize me */
                        #ifdef MPI_PARALLEL
                        ambient_assert(false);
                        ambient::push(ambient::reshape_r2l_l_kernel, ambient::reshape_r2l_c_kernel, in_block, out_block, 
                                      out_left_offset, in_right_offset, physical_i[s].second, left_i[l].second, right_i[r].second);
                        #else
                        for (size_t ss = 0; ss < physical_i[s].second; ++ss)
                            for (size_t rr = 0; rr < right_i[r].second; ++rr)
                                memcpy(&out_block(out_left_offset + ss*left_i[l].second, rr),
                                       &in_block(0, in_right_offset + ss*right_i[r].second+rr),
                                       sizeof(typename Matrix::value_type) * left_i[l].second);
                                /*for(size_t ll = 0; ll < left_i[l].second; ++ll)
                                      out_block(out_left_offset + ss*left_i[l].second+ll, rr) = 
                                      in_block(ll, in_right_offset + ss*right_i[r].second+rr);*/
                        #endif
                    }

                    if (pretend)
                        m2.reserve(out_l_charge, out_r_charge,
                                   out_left_offset + physical_i[s].second * left_i[l].second, right_i[r].second);
                }
            }
        }
    }
    
    timer.end();
    
//    assert(m2.right_basis() == right_i);
}

#endif /* RESHAPE_H */

