#ifndef RESHAPE_H
#define RESHAPE_H

#include "block_matrix/indexing.h"

template<class SymmGroup, unsigned long L>
std::pair<typename SymmGroup::charge, std::size_t>
calculate_index(boost::array<Index<SymmGroup>, L> const & dims,
                boost::array<std::pair<typename SymmGroup::charge, std::size_t>, L> const & idx)
{
    typedef typename SymmGroup::charge charge;
    using std::size_t;
    
    std::pair<charge, size_t> ret = std::make_pair(SymmGroup::SingletCharge, 0);
    
    // first step: calculate final charge
    for (size_t k = 0; k < L; ++k)
        ret.first = SymmGroup::fuse(ret.first, idx[k].first);
    
    // second step: strides
    boost::array<size_t, L> strides;
    strides[L-1] = 1;
    for (size_t k = L-1; k > 0; --k)
        strides[k-1] = strides[k] * dims[k].size_of_block(idx[k].first);
    
//    std::copy(strides.begin(), strides.end(), std::ostream_iterator<size_t>(cout, " ")); cout << endl;
    
    // last step: index
    for (size_t k = 0; k < L; ++k)
        ret.second += strides[k] * idx[k].second;
    
    return ret;
}

template<class Matrix, class SymmGroup>
void reshape_left_to_right(Index<SymmGroup> physical_i,
                           Index<SymmGroup> left_i,
                           Index<SymmGroup> right_i,
                           block_matrix<Matrix, SymmGroup> const & m1,
                           block_matrix<Matrix, SymmGroup> & m2)
{
    using std::size_t;
    
//    assert(m1.left_basis() == physical_i*left_i);
//    assert(m1.right_basis() == right_i);
    
    Index<SymmGroup> new_phys_i = adjoin(physical_i);
    Index<SymmGroup> newl = left_i, newr = new_phys_i * right_i;
    common_subset(newl, newr);
    m2 = block_matrix<Matrix, SymmGroup>(newl, newr);
    
//    cout << "Reshaping l2r: " << endl;
//    cout << physical_i << " " << left_i << " " << right_i << endl;
//    cout << physical_i * left_i << " " << right_i << endl;
//    cout << newl << " " << newr << endl;
    
    typedef typename Index<SymmGroup>::basis_iterator bit;
    
    for (bit s = physical_i.basis_begin(); !s.end(); ++s)
        for (bit l = left_i.basis_begin(); !l.end(); ++l) {
            std::pair<typename SymmGroup::charge, std::size_t>
            i1 = calculate_index(physical_i ^ left_i, *s ^ *l);
            
            for (bit r = right_i.basis_begin(); !r.end(); ++r) {
                if (i1.first != (*r).first)
                    continue;
                
//                cout << "( " << i1.first << "," << (*r).first << " ) -> ( ";
//                cout << (*l).first << "," << 
//                calculate_index(new_phys_i ^ right_i, std::make_pair(-(*s).first, (*s).second) ^ *r).first << ")" << endl;
                m2(*l, calculate_index(new_phys_i ^ right_i, std::make_pair(-(*s).first, (*s).second) ^ *r)) = m1(i1, *r);
            }
        }
}

template<class Matrix, class SymmGroup>
void reshape_right_to_left(Index<SymmGroup> physical_i,
                           Index<SymmGroup> left_i,
                           Index<SymmGroup> right_i,
                           block_matrix<Matrix, SymmGroup> const & m1,
                           block_matrix<Matrix, SymmGroup> & m2)
{
    using std::size_t;
    
//    assert(m1.left_basis() == left_i);
//    assert(m1.right_basis() == physical_i*right_i);
    
    Index<SymmGroup> new_phys_i = adjoin(physical_i);
    Index<SymmGroup> newl = new_phys_i * left_i, newr = right_i;
    common_subset(newl, newr);
    m2 = block_matrix<Matrix, SymmGroup>(newl, newr);
    
    typedef typename Index<SymmGroup>::basis_iterator bit;
    
    for (bit s = physical_i.basis_begin(); !s.end(); ++s)
        for (bit l = left_i.basis_begin(); !l.end(); ++l) {
            std::pair<typename SymmGroup::charge, std::size_t>
            i1 = calculate_index(physical_i ^ left_i, *s ^ *l);
            
            for (bit r = right_i.basis_begin(); !r.end(); ++r) {
                if (i1.first != (*r).first)
                    continue;
                
                m2(i1, *r) = m1(*l, calculate_index(physical_i ^ right_i, std::make_pair(-(*s).first, (*s).second) ^ *r));
            }
        }
}

#endif /* RESHAPE_H */
