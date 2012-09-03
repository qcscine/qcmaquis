/*
*Very Large Integer Library, License - Version 1.0 - May 3rd, 2012
*
*Timothee Ewart - University of Geneva, 
*Andreas Hehn - Swiss Federal Institute of technology Zurich.
*
*Permission is hereby granted, free of charge, to any person or organization
*obtaining a copy of the software and accompanying documentation covered by
*this license (the "Software") to use, reproduce, display, distribute,
*execute, and transmit the Software, and to prepare derivative works of the
*Software, and to permit third-parties to whom the Software is furnished to
*do so, all subject to the following:
*
*The copyright notices in the Software and this entire statement, including
*the above license grant, this restriction and the following disclaimer,
*must be included in all copies of the Software, in whole or in part, and
*all derivative works of the Software, unless such copies or derivative
*works are solely in the form of machine-executable object code generated by
*a source language processor.
*
*THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
*SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
*FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
*ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
*DEALINGS IN THE SOFTWARE.
*/

#include "vli/detail/cpu/x86_64/kernel_implementation_macros.h"
#include <limits>

namespace vli {
    namespace detail {
                    // new functions type : VLI<n*64> - VLI<n*64> : sub128_128, sub192_192 ...
                    template <std::size_t NumWords>
                    void sub(boost::uint64_t* x,  boost::uint64_t const*  y);
                     
                     #define FUNCTION_sub_nbits_nbits(z, n, unused) \
                         template<> \
                         void sub<(n+2)>( boost::uint64_t* x,  boost::uint64_t const* y){\
                             boost::uint32_t const* cy = const_cast<boost::uint32_t*>(reinterpret_cast<boost::uint32_t const*>(y)); \
                             boost::uint32_t* cx = reinterpret_cast<boost::uint32_t*>(x);                                           \
                             boost::uint32_t borrow(0);                                                                             \
                             for(int i(0); i<2*(n+2);++i){                                                                          \
                                boost::uint64_t tmp = static_cast<uint64_t>(cx[i]) - cy[i] - borrow;                                \
                                borrow = (tmp >> (std::numeric_limits<boost::uint64_t>::digits)-1) ;                                \
                                cx[i] = tmp;                                                                                        \
                             }                                                                                                      \
                          }                                                                                                         \

                     BOOST_PP_REPEAT(VLI_MAX_ITERATION, FUNCTION_sub_nbits_nbits, ~)
                     #undef FUNCTION_sub_nbits_nbits
                     
                    template <std::size_t NumWords>
                    void sub(boost::uint64_t* x,  boost::uint64_t const b);

                     //new functions type : VLI<n*64> - VLI<64> : sub192_64, sub256_64
                     //the case is done after sub128_64
                     #define FUNCTION_sub_nbits_64bits(z, n, unused) \
                          template<> \
                          void sub<(n+2)>( boost::uint64_t* x,  boost::uint64_t const b){    \
                            boost::uint32_t cb = static_cast<boost::uint32_t>(b);                               \
                            boost::uint32_t* cx = reinterpret_cast<boost::uint32_t*>(x);                        \
                            boost::uint32_t sign = cb >> std::numeric_limits<boost::uint32_t>::digits - 1;      \
                            sign = -sign;                                                                       \
                            boost::uint32_t borrow(0);                                                          \
                            boost::uint64_t tmp = static_cast<uint64_t>(cx[0]) - cb;                            \
                            cx[0] = tmp;                                                                        \
                            borrow = tmp >> (std::numeric_limits<boost::uint64_t>::digits-1);                   \
                            for(int i(1); i<2*(n+2);++i){                                                       \
                               boost::uint64_t tmp = static_cast<uint64_t>(cx[i]) - sign - borrow;              \
                               borrow = tmp >> (std::numeric_limits<boost::uint64_t>::digits-1);                \
                               cx[i] = tmp;                                                                     \
                            }                                                                                   \
                          }                                                                                     \

                     BOOST_PP_REPEAT(VLI_MAX_ITERATION, FUNCTION_sub_nbits_64bits, ~)
                     #undef FUNCTION_sub_nbits_64bits

    } // end namespace detail
} // end namespace vli
