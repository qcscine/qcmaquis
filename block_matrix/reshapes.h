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
    
    m2 = block_matrix<Matrix, NullGroup>(physical_i*left_i, right_i);
    
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
    
    m2 = block_matrix<Matrix, NullGroup>(left_i, physical_i*right_i);
    
    for (size_t s = 0; s < physical_i[0].second; ++s)
        for (size_t l = 0; l < left_i[0].second; ++l)
            for (size_t r = 0; r < right_i[0].second; ++r)
                m2[0](s*left_i[0].second+l, r) = m1[0](l, s*right_i[0].second+r);
}
                           

#endif /* RESHAPE_H */
