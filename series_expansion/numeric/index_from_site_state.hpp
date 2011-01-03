#ifndef HP2C__INDEX_FROM_SITESTATE_HPP
#define HP2C__INDEX_FROM_SITESTATE_HPP
namespace series_expansion
{


//Helper Metaprogram
namespace detail
{
    template <unsigned int StatesPerSite>
    struct bits_per_site_helper
    {
        enum { value = 1+ bits_per_site_helper< (StatesPerSite>>1) >::value };
    };

    template <>
    struct bits_per_site_helper<0>
    {
        enum { value = 0 };
    };
    
    template <unsigned int StatesPerSite>
    struct bits_per_site
    {
        enum { value = bits_per_site_helper<StatesPerSite-1>::value };
    };
}

template <unsigned int NumSiteStates, typename Index>
const Index index_from_site_state_pair(Index const& base_index, unsigned int site1, unsigned int site_state1, unsigned int site2=0, unsigned int site_state2=0)
{
    assert(site_state1 < NumSiteStates);
    assert(site_state2 < NumSiteStates);
    // Generate a bit mask for 4 site states eg. 0b...00100 - 1 = 0b...000011 
    const Index mask((1<<detail::bits_per_site<NumSiteStates>::value) - 1);
    Index i = (~(mask<<(site2 * detail::bits_per_site<NumSiteStates>::value )) & base_index) | site_state2<<(site2*detail::bits_per_site<NumSiteStates>::value);
    return (~(mask<<(site1 * detail::bits_per_site<NumSiteStates>::value )) & i) | site_state1<<(site1*detail::bits_per_site<NumSiteStates>::value);
}

template <unsigned int NumSiteStates, typename Index>
unsigned int get_site_state_from_index(Index const& i, unsigned int site)
{
    const Index mask((1<<detail::bits_per_site<NumSiteStates>::value) - 1);
    return mask & static_cast<Index>(i>>(site*detail::bits_per_site<NumSiteStates>::value));
}
}


#endif //HP2C__INDEX_FROM_SITESTATE_HPP
