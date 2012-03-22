#ifndef VLI_BIT_MASKS_HPP
#define VLI_BIT_MASKS_HPP

#include <boost/static_assert.hpp>

namespace vli
{

// LOG_BASE
template <typename T>
struct data_bits
{
    //The following line ensures that the default template is never used
    // -> only specializations are valid.
    BOOST_STATIC_ASSERT(sizeof(T) == 0); 
    enum { value = 0 };
};


template <>
struct data_bits<unsigned long int>
{
   BOOST_STATIC_ASSERT( sizeof(long unsigned int) == 8 ); 
    enum { value = 63 };
};

template <>
struct data_bits<unsigned int>
{
    BOOST_STATIC_ASSERT( sizeof(unsigned int) == 4 );
    enum { value = 30 };
};

template <typename T>
struct base
{
    enum { value = static_cast<T>(1)<<(data_bits<T>::value)};
};

template <typename T>
struct data_mask
{
    enum { value = static_cast<T>(base<T>::value)-1};
};

template <typename T>
struct max_value
{
    enum { value = data_mask<T>::value };
};

} //namespace vli

#endif //VLI_BIT_MASKS_HPP
