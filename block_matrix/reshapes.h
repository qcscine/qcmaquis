#ifndef RESHAPE_H
#define RESHAPE_H

// A( [sigma, l] [r] ) -> A( [l] [sigma, r] )

template<class Matrix>
void reshape_left_to_right(Index<NullGroup> physical_i,
                           Index<NullGroup> left_i,
                           Index<NullGroup> right_i,
                           block_matrix<Matrix, NullGroup> const & m1,
                           block_matrix<Matrix, NullGroup> & m2)
{
    using std::size_t;
    
    assert(m1.left_basis() == physical_i*left_i);
    assert(m1.right_basis() == right_i);
    
    m2 = block_matrix<Matrix, NullGroup>(left_i, physical_i*right_i);
    
    for (size_t s = 0; s < physical_i[0].second; ++s)
        for (size_t l = 0; l < left_i[0].second; ++l)
            for (size_t r = 0; r < right_i[0].second; ++r)
                m2[0](l, s*right_i[0].second+r) = m1[0](s*left_i[0].second+l, r);
}

template<class Matrix>
void reshape_right_to_left(Index<NullGroup> physical_i,
                           Index<NullGroup> left_i,
                           Index<NullGroup> right_i,
                           block_matrix<Matrix, NullGroup> const & m1,
                           block_matrix<Matrix, NullGroup> & m2)
{
    using std::size_t;
    
    assert(m1.left_basis() == left_i);
    assert(m1.right_basis() == physical_i*right_i);
    
    m2 = block_matrix<Matrix, NullGroup>(physical_i*left_i, right_i);
    
    for (size_t s = 0; s < physical_i[0].second; ++s)
        for (size_t l = 0; l < left_i[0].second; ++l)
            for (size_t r = 0; r < right_i[0].second; ++r)
                m2[0](s*left_i[0].second+l, r) = m1[0](l, s*right_i[0].second+r);
}
                           
template<class SymmGroup, int L>
std::pair<typename SymmGroup::charge, std::size_t>
calculate_index(boost::array<Index<SymmGroup>, L> const & dims,
                boost::array<std::pair<typename SymmGroup::charge, std::size_t>, L> const & idx)
{
    typedef typename SymmGroup::charge charge;
    using std::size_t;
    
    std::pair<charge, size_t> ret = make_pair(SymmGroup::SingletCharge, 0);
    
    // first step: calculate final charge
    for (size_t k = 0; k < L; ++k)
        ret.first = SymmGroup::fuse(ret.first, idx[k].first);
    
    // second step: strides
    boost::array<size_t, L> strides;
    for (size_t k = 0; k < L; ++k) {
        if (k == 0)
            strides[0] = dims[k].get_size(idx[k].first);
        else
            strides[k] = strides[k-1] * dims[k].get_size(idx[k].first);
    }
    
    // last step: index
    for (size_t k = 0; k < L; ++k)
        ret.second += strides[k] * idx[k].second;
    
    return ret;
}

template<class SymmGroup>
boost::array<Index<SymmGroup>, 2> operator,(Index<SymmGroup> const & i1,
                                            Index<SymmGroup> const & i2)
{
    boost::array<Index<SymmGroup>, 2> ret;
    ret[0] = i1;
    ret[1] = i2;
    return ret;
}

template<class SymmGroup, int L>
boost::array<Index<SymmGroup>, L+1> operator,(boost::array<Index<SymmGroup>, L> const & i1,
                                              Index<SymmGroup> const & i2)
{
    boost::array<Index<SymmGroup>, L+1> ret;
    std::copy(i1.begin(), i1.end(), ret.begin());
    ret[L] = i2;
    return ret;
}

template<class SymmGroup>
boost::array<std::pair<typename SymmGroup::charge, std::size_t>, 2>
operator,(std::pair<typename SymmGroup::charge, std::size_t> const & i1,
          std::pair<typename SymmGroup::charge, std::size_t> const & i2)
{
    boost::array<std::pair<typename SymmGroup::charge, std::size_t>, 2> ret;
    ret[0] = i1;
    ret[1] = i2;
    return ret;
}

template<class SymmGroup, int L>
boost::array<std::pair<typename SymmGroup::charge, std::size_t>, L+1>
operator,(boost::array<std::pair<typename SymmGroup::charge, std::size_t>, L> const & i1,
          std::pair<typename SymmGroup::charge, std::size_t> const & i2)
{
    boost::array<std::pair<typename SymmGroup::charge, std::size_t>, L+1> ret;
    std::copy(i1.begin(), i1.end(), ret.begin());
    ret[L] = i2;
    return ret;
}

#endif /* RESHAPE_H */
