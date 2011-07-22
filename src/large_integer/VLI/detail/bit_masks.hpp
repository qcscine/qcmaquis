#ifndef VLI_BIT_MASKS_HPP
#define VLI_BIT_MASKS_HPP

#include <boost/static_assert.hpp>

namespace vli
{

// LOG_BASE
template <typename T>
struct data_bits
{
//    BOOST_STATIC_ASSERT( 1 == 0);
//    BOOST_STATIC_ASSERT_MSG( false, "VLI can not use this type as BaseInt!");
    enum { value = 0 };
};

// BASE
template <typename T>
struct base
{
    enum { value = static_cast<T>(1)<<data_bits<T>::value };
};

// BASE_HALF
template <typename T>
struct base_half
{
    BOOST_STATIC_ASSERT( data_bits<T>::value % 2 == 0 );
    enum { value = static_cast<T>(1)<<(data_bits<T>::value/2) };
};

// BASE_MINUS
template <typename T>
struct data_mask
{
    enum { value = static_cast<T>(base<T>::value)-1 };
};

// MASK_DOWN
template <typename T>
struct mask_down
{
    enum { value = static_cast<T>(base_half<T>::value)-1 };
};

// MASK_UP
template <typename T>
struct mask_up
{
    enum { value = (static_cast<T>(data_mask<T>::value) ^ static_cast<T>(mask_down<T>::value)) };
};

template <typename T>
struct max_value
{
    enum { value = data_mask<T>::value };
};



template <>
struct data_bits<unsigned long int>
{
   BOOST_STATIC_ASSERT( sizeof(long unsigned int) == 8 ); 
    enum { value = 62 };
};

template <>
struct data_bits<unsigned int>
{
    BOOST_STATIC_ASSERT( sizeof(unsigned int) == 4 );
    enum { value= 30 };
};

template <>
struct data_bits<int>
{
    BOOST_STATIC_ASSERT( sizeof(int) == 4 );
    enum { value= 30 };
};

} //namespace vli

#endif //VLI_BIT_MASKS_HPP
