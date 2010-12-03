#ifndef RESHAPE_H
#define RESHAPE_H

#include "indexing.h"

template<class SymmGroup, int L>
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
        strides[k-1] = strides[k] * dims[k].get_size(idx[k].first);
    
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
    
    assert(m1.left_basis() == physical_i*left_i);
    assert(m1.right_basis() == right_i);
    
    m2 = block_matrix<Matrix, NullGroup>(left_i, physical_i*right_i);
    
    typedef typename Index<SymmGroup>::basis_iterator bit;
    
    for (bit s = physical_i.basis_begin(); !s.end(); ++s)
        for (bit l = left_i.basis_begin(); !l.end(); ++l)
            for (bit r = right_i.basis_begin(); !r.end(); ++r)
                m2(calculate_index<SymmGroup, 1>(_(left_i), _(*l)),
                   calculate_index<SymmGroup, 2>(physical_i ^ right_i,
                                                 *s ^ *r))
                =
                m1(calculate_index<SymmGroup, 2>(physical_i ^ left_i,
                                                 *s ^ *l),
                   calculate_index<SymmGroup, 1>(_(right_i), _(*r)));
}

template<class Matrix, class SymmGroup>
void reshape_right_to_left(Index<SymmGroup> physical_i,
                           Index<SymmGroup> left_i,
                           Index<SymmGroup> right_i,
                           block_matrix<Matrix, SymmGroup> const & m1,
                           block_matrix<Matrix, SymmGroup> & m2)
{
    using std::size_t;
    
    assert(m1.left_basis() == left_i);
    assert(m1.right_basis() == physical_i*right_i);
    
    m2 = block_matrix<Matrix, NullGroup>(physical_i*left_i, right_i);
    
    typedef typename Index<SymmGroup>::basis_iterator bit;
    
    for (bit s = physical_i.basis_begin(); !s.end(); ++s)
        for (bit l = left_i.basis_begin(); !l.end(); ++l)
            for (bit r = right_i.basis_begin(); !r.end(); ++r)
                m2(calculate_index<SymmGroup, 2>(physical_i ^ left_i,
                                                 *s ^ *l),
                   calculate_index<SymmGroup, 1>(_(right_i), _(*r)))
                =
                m1(calculate_index<SymmGroup, 1>(_(left_i), _(*l)),
                   calculate_index<SymmGroup, 2>(physical_i ^ right_i,
                                                 *s ^ *r));
}

#endif /* RESHAPE_H */
