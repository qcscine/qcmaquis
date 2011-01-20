#ifndef RESHAPE_H
#define RESHAPE_H

#include <map>

#include "block_matrix/indexing.h"

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
        std::map<charge, size_t> in_offsets, out_offsets;
        
        if (run == 1)
            m2.allocate_blocks();
        
        for (size_t s = 0; s < physical_i.size(); ++s)
            for (size_t l = 0; l < left_i.size(); ++l)
                for (size_t r = 0; r < right_i.size(); ++r)
                {
                    bool pretend = (run == 0);
                    
                    charge in_l_charge = SymmGroup::fuse(physical_i[s].first, left_i[l].first);
                    charge in_r_charge = right_i[r].first;
                    charge out_l_charge = left_i[l].first;
                    charge out_r_charge = SymmGroup::fuse(-physical_i[s].first, right_i[r].first);

                    if (! m1.has_block(in_l_charge, in_r_charge) )
                        continue;
                    
                    size_t in_left_offset = in_offsets[in_l_charge];
                    size_t out_right_offset = out_offsets[out_r_charge];
                    
                    if (!pretend) {
                        Matrix const & in_block = m1(in_l_charge, in_r_charge);
                        Matrix & out_block = m2(out_l_charge, out_r_charge);
                        
                        /* optimize me */
                        for (size_t ss = 0; ss < physical_i[s].second; ++ss)
                            for (size_t rr = 0; rr < right_i[r].second; ++rr)
                                for (size_t ll = 0; ll < left_i[l].second; ++ll)
                                    out_block(ll, out_right_offset + ss*right_i[r].second+rr) = in_block(in_left_offset + ss*left_i[l].second+ll, rr);
                    }
                        
                    in_offsets[in_l_charge] += left_i[l].second * physical_i[s].second;
                    out_offsets[out_r_charge] += right_i[r].second * physical_i[s].second;

                    if (pretend)
                        m2.reserve(out_l_charge, out_r_charge,
                                   left_i[l].second, out_right_offset + physical_i[s].second * right_i[r].second);
                }
    }
    
    timer.end();
    
//    assert(m2.left_basis() == left_i);
}

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
        std::map<charge, size_t> in_offsets, out_offsets;
        
        if (run == 1)
            m2.allocate_blocks();
    
        for (size_t s = 0; s < physical_i.size(); ++s)
            for (size_t l = 0; l < left_i.size(); ++l)
                for (size_t r = 0; r < right_i.size(); ++r)
                {
                    bool pretend = (run == 0);
                    
                    charge in_l_charge = left_i[l].first;
                    charge in_r_charge = SymmGroup::fuse(-physical_i[s].first, right_i[r].first);
                    charge out_l_charge = SymmGroup::fuse(physical_i[s].first, left_i[l].first);
                    charge out_r_charge = right_i[r].first;
                    
                    if (! m1.has_block(in_l_charge, in_r_charge) )
                        continue;
                    
                    size_t in_right_offset = in_offsets[in_r_charge];
                    size_t out_left_offset = out_offsets[out_l_charge];
                    
                    if (!pretend) {
                        Matrix const & in_block = m1(in_l_charge, in_r_charge);
                        Matrix & out_block = m2(out_l_charge, out_r_charge);
                        
                        /* optimize me */
                        for (size_t ss = 0; ss < physical_i[s].second; ++ss)
                            for (size_t rr = 0; rr < right_i[r].second; ++rr)
                                memcpy(&out_block(out_left_offset + ss*left_i[l].second, rr),
                                       &in_block(0, in_right_offset + ss*right_i[r].second+rr),
                                       sizeof(typename Matrix::value_type) * left_i[l].second);
                                /*
                                for (size_t ll = 0; ll < left_i[l].second; ++ll)
                                    out_block(out_left_offset + ss*left_i[l].second+ll, rr) 
                                 = in_block(ll, in_right_offset + ss*right_i[r].second+rr);*/
                    }
                        
                    in_offsets[in_r_charge] += right_i[r].second * physical_i[s].second;
                    out_offsets[out_l_charge] += left_i[l].second * physical_i[s].second;

                    if (pretend)
                        m2.reserve(out_l_charge, out_r_charge,
                                   out_left_offset + physical_i[s].second * left_i[l].second, right_i[r].second);
                }
    }
    
    timer.end();
    
//    assert(m2.right_basis() == right_i);
}

#endif /* RESHAPE_H */

