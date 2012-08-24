/*
*Very Large Integer Library, License - Version 1.0 - May 3rd, 2012
*
*Timothee Ewart - University of Geneva, 
*Andreas Hehn - Swiss Federal Institute of technology Zurich.
*Maxim Milakov – NVIDIA
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

 #ifndef VLI_VARIABLES_GPU_H
 #define VLI_VARIABLES_GPU_H

namespace vli {

    namespace detail {


    template <int Var, int NumVars, unsigned int Order>
    struct result_stride {
        static unsigned int const value = stride<Var,NumVars,2*Order>::value;
    };

    template <int Var, int NumVars, unsigned int Order>
    struct stride_pad {
        static unsigned int const value = Var < NumVars ? Order+1 : 0;
    };


    struct SumBlockSize {
       enum { value = 256};
    };

//    TODO Remove this, all num_coefficients in the gpu code should probably be replaced by num_coefficients<...<2*Order>,...>
//
//    template<class MaxOrder, int NumVars>
//    struct num_coefficients;
//
//        
//    template<class MaxOrder>
//    struct num_coefficients<MaxOrder, 0>{
//        enum {value = 1};
//    };
//
//    template<unsigned int Order, int NumVars>
//    struct num_coefficients<max_order_each<Order>, NumVars>{
//        enum {value = (2*Order+1)*num_coefficients<max_order_each<Order>, NumVars-1>::value};
//    };
//
//    template<unsigned int Order, int NumVars>
//    struct num_coefficients<max_order_combined<Order>, NumVars>{
//        enum {value = vli::detail::max_order_combined_helpers::size<NumVars+1, 2*Order>::value};
//    };

    template<class MaxOrder, int NumVars>
    struct mul_block_size {
        enum {value = num_coefficients<MaxOrder,NumVars>::value/2 >= 256U) ? 256U
               : (num_coefficients<MaxOrder,NumVars>::value/2U+32U-1U)/32U*32U };
    };

 
    template<class MaxOrder, int NumVars>
    struct MaxIterationCount {
        enum {value = (num_coefficients<MaxOrder,NumVars>::value+32U-1U)/32U};
    };

    // replace in the code MaxNumberCoefficientExtend by num_coefficients

    
    template<std::size_t Size>
    struct size_pad{
        enum {value = Size | 1};
    };

    }
}

#endif
