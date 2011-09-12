#ifndef VLI_NUMBER_TRAITS_HPP
#define VLI_NUMBER_TRAITS_HPP

#include "detail/bit_masks.hpp"

namespace vli
{

template <typename Vli>
struct max_int_value
{
    enum { value = data_mask<typename Vli::value_type>::value };
};

} // namespace vli

#endif //VLI_NUMBER_TRAITS_HPP
