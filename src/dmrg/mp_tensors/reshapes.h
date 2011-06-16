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
                        ambient::push(ambient::reshape_l2r_l, ambient::reshape_l2r_c, in_block, out_block, 
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
                        ambient::push(ambient::reshape_r2l_l, ambient::reshape_r2l_c, out_block, in_block, 
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


template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> reshape_2site_op (Index<SymmGroup> const & phys,
                                                  block_matrix<Matrix, SymmGroup> const & A)
{
    typedef typename SymmGroup::charge charge;
    block_matrix<Matrix, SymmGroup> ret;
    
    ProductBasis<SymmGroup> pb(phys, phys);
    ProductBasis<SymmGroup> pb_new(phys, adjoin(phys));
    
    typedef typename Index<SymmGroup>::basis_iterator bi_t;
    for (bi_t s1 = phys.basis_begin(); !s1.end(); ++s1)
        for (bi_t s2 = phys.basis_begin(); !s2.end(); ++s2)
            for (bi_t s3 = phys.basis_begin(); !s3.end(); ++s3)
                for (bi_t s4 = phys.basis_begin(); !s4.end(); ++s4)
                {                    
                    charge in_left_c = SymmGroup::fuse(s1->first, s2->first);
                    charge in_right_c = SymmGroup::fuse(s3->first, s4->first);
                    
                    if (!A.has_block(in_left_c, in_right_c))
                        continue;
                    
                    charge out_left_c = SymmGroup::fuse(s1->first, -s3->first);
                    charge out_right_c = SymmGroup::fuse(s4->first, -s2->first);
                    
                    if (!ret.has_block(out_left_c, out_right_c))
                        ret.insert_block(Matrix(pb_new.size(s1->first, -s3->first),
                                                pb_new.size(s4->first, -s2->first),
                                                0),
                                         out_left_c, out_right_c);
                    
                    std::size_t in_left_offset = pb(s1->first, s2->first);
                    std::size_t in_right_offset = pb(s3->first, s4->first);
                    
                    std::size_t out_left_offset = pb_new(s1->first, -s3->first);
                    std::size_t out_right_offset = pb_new(s4->first, -s2->first);
                    
                    ret(std::make_pair(out_left_c, out_left_offset + s1->second*phys.size_of_block(s3->first)+s3->second),
                        std::make_pair(out_right_c, out_right_offset + s2->second*phys.size_of_block(s4->first)+s4->second))
                    = A(std::make_pair(in_left_c, in_left_offset + s1->second*phys.size_of_block(s2->first) + s2->second),
                        std::make_pair(in_right_c, in_right_offset + s3->second*phys.size_of_block(s4->first) + s4->second));
                    
                }
    
    // Removing empty blocks
    Index<U1> out_left = phys*adjoin(phys);
    Index<U1> out_right = phys*adjoin(phys);
    for (typename Index<U1>::const_iterator it1 = out_left.begin(); it1 != out_left.end(); it1++)
    {
        for (typename Index<U1>::const_iterator it2 = out_right.begin(); it2 != out_right.end(); it2++)
        {
            bool empty = true;
            
            if (!ret.has_block(it1->first, it2->first)) continue;
            
            Matrix tmp = ret(it1->first, it2->first);
            for (int i=0; i<blas::num_rows(tmp); ++i)
            {
                for (int j=0; j<blas::num_cols(tmp); ++j)
                    if (tmp(i,j) != 0) {
                        empty=false;
                        break;
                    }
                if (!empty) break;
            }
            
            if (empty)
                ret.remove_block(it1->first, it2->first);
        }
    }
    
    return ret;
}

template<class Matrix, class SymmGroup>
std::vector<block_matrix<Matrix, SymmGroup> > reshape_right_to_list (Index<SymmGroup> const & phys,
                                                                     block_matrix<Matrix, SymmGroup> const & A)
{
    typedef typename SymmGroup::charge charge;
    std::vector<block_matrix<Matrix, SymmGroup> > ret;
    
    Index<SymmGroup> aux_i = A.right_basis();
    ProductBasis<SymmGroup> pb(phys, adjoin(phys));

    typedef typename Index<SymmGroup>::basis_iterator bi_t;
    for (bi_t b = aux_i.basis_begin(); !b.end(); ++b)
    {
    	block_matrix<Matrix, SymmGroup> Ai;
        for (bi_t s1 = phys.basis_begin(); !s1.end(); ++s1)
            for (bi_t s2 = phys.basis_begin(); !s2.end(); ++s2)
            {                    
                charge in_left_c = SymmGroup::fuse(s1->first, -s2->first);
                charge in_right_c = b->first;
                
                if (!A.has_block(in_left_c, in_right_c))
                    continue;
                
                charge out_left_c = s1->first;
                charge out_right_c = s2->first;
                
                if (!Ai.has_block(out_left_c, out_right_c))
                    Ai.insert_block(Matrix(phys.size_of_block(s1->first),
                                            phys.size_of_block(s2->first),
                                            0),
                                     out_left_c, out_right_c);
                
                std::size_t in_left_offset = pb(s1->first, -s2->first);
                std::size_t in_right_offset = 0;
                
                std::size_t out_left_offset = 0;
                std::size_t out_right_offset = 0;
                
                Ai(*s1, *s2)
                = A(std::make_pair(in_left_c, in_left_offset + s1->second*phys.size_of_block(s2->first) + s2->second),
                    *b);
                
            }

        // Removing empty blocks
        for (int n=0; n<Ai.n_blocks(); ++n)
        {
        	bool empty = true;
        	for (int i=0; i<blas::num_rows(Ai[n]) && empty; ++i)
        		for (int j=0; j<blas::num_cols(Ai[n]) && empty; ++j)
        			if (Ai[n](i,j) != 0)
        				empty=false;

        	if (empty)
        		Ai.remove_block(n);
        }

        ret.push_back(Ai);
    }
    assert( ret.size() == aux_i.sum_of_sizes() );

    return ret;
}

template<class Matrix, class SymmGroup>
std::vector<block_matrix<Matrix, SymmGroup> > reshape_left_to_list (Index<SymmGroup> const & phys,
                                                                    block_matrix<Matrix, SymmGroup> const & A)
{
	typedef typename SymmGroup::charge charge;
	std::vector<block_matrix<Matrix, SymmGroup> > ret;

	Index<SymmGroup> aux_i = A.left_basis();
	ProductBasis<SymmGroup> pb(phys, adjoin(phys));

	typedef typename Index<SymmGroup>::basis_iterator bi_t;
	for (bi_t b = aux_i.basis_begin(); !b.end(); ++b)
	{
		block_matrix<Matrix, SymmGroup> Ai;
		for (bi_t s1 = phys.basis_begin(); !s1.end(); ++s1)
			for (bi_t s2 = phys.basis_begin(); !s2.end(); ++s2)
			{
				charge in_right_c = SymmGroup::fuse(s2->first, -s1->first);
				charge in_left_c = b->first;

				if (!A.has_block(in_left_c, in_right_c))
					continue;

				charge out_left_c = s1->first;
				charge out_right_c = s2->first;

				if (!Ai.has_block(out_left_c, out_right_c))
					Ai.insert_block(Matrix(phys.size_of_block(s1->first),
							phys.size_of_block(s2->first),
							0),
							out_left_c, out_right_c);

				std::size_t in_left_offset = 0;
				std::size_t in_right_offset = pb(s2->first, -s1->first);

				std::size_t out_left_offset = 0;
				std::size_t out_right_offset = 0;

				Ai(std::make_pair(out_left_c, out_left_offset + s1->second),
						std::make_pair(out_right_c, out_right_offset + s2->second))
				= A(std::make_pair(in_left_c, in_left_offset + b->second),
						std::make_pair(in_right_c, in_right_offset + s1->second*phys.size_of_block(s2->first) + s2->second));

			}

		// Removing empty blocks
		for (int n=0; n<Ai.n_blocks(); ++n)
		{
			bool empty = true;
			for (int i=0; i<blas::num_rows(Ai[n]) && empty; ++i)
				for (int j=0; j<blas::num_cols(Ai[n]) && empty; ++j)
					if (Ai[n](i,j) != 0)
						empty=false;

			if (empty)
				Ai.remove_block(n);
		}

		ret.push_back(Ai);
	}
	assert( ret.size() == aux_i.sum_of_sizes() );

	return ret;
}


#endif /* RESHAPE_H */

