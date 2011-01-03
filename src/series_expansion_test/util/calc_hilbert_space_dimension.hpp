#ifndef HP2C__CALC_HILBERT_SPACE_DIMENSION_HPP
#define HP2C__CALC_HILBERT_SPACE_DIMENSION_HPP

namespace series_expansion
{

template <typename SizeType, typename Operator, typename Graph>
SizeType calc_hilbert_space_dimension(SizeType const& s, Operator const& Ho, Graph const& g)
{
    SizeType dim(1);
    unsigned int num_vtcs = num_vertices(g);
    for(unsigned int i=0; i < num_vtcs; ++i)
    {
        dim *= Operator::states_per_site;
    }
    return dim;
}

}
#endif //HP2C__CALC_HIBLERT_SPACE_DIMENSION_HPP
