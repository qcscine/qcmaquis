/*****************************************************************************
 *
 * MAQUIS MATRICES Project
 *
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 * 
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 * 
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not, the license is also
 * available from http://alps.comp-phys.org/.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
 * DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

#ifndef MAQUIS_FUNCTION_OBJECTS_H
#define MAQUIS_FUNCTION_OBJECTS_H

#include <alps/numeric/conj.hpp>

namespace utils {
    using alps::numeric::conj;
    
#define DEFINE_FUNCTION_OBJECT(name, return_type, arg_type) \
struct functor_##name { template<class T> return_type operator() (arg_type t) { return name(t); } };

#define DEFINE_VOID_FUNCTION_OBJECT(name, arg_type) \
struct functor_##name { template<class T> void operator() (arg_type t) { name(t); } };
    
    DEFINE_FUNCTION_OBJECT(trace, typename maquis::traits::scalar_type<T>::type, T const &)
    DEFINE_FUNCTION_OBJECT(norm_square, typename maquis::traits::real_type<T>::type, T const &)
    DEFINE_FUNCTION_OBJECT(transpose, T, T const &)
    DEFINE_VOID_FUNCTION_OBJECT(transpose_inplace, T &)
    DEFINE_FUNCTION_OBJECT(conj, T, T)
    DEFINE_VOID_FUNCTION_OBJECT(conj_inplace, T&)
    DEFINE_FUNCTION_OBJECT(conjugate, T, T const &)
    DEFINE_FUNCTION_OBJECT(adjoint, T, T const &)
    DEFINE_VOID_FUNCTION_OBJECT(adjoint_inplace, T &)
    DEFINE_FUNCTION_OBJECT(sqrt, T, T const &)
    DEFINE_VOID_FUNCTION_OBJECT(sqrt_inplace, T &)
    DEFINE_FUNCTION_OBJECT(exp, T, T const &)
    DEFINE_FUNCTION_OBJECT(size_of, std::size_t, T const &)
    
#undef DEFINE_FUNCTION_OBJECT
#undef DEFINE_VOID_FUNCTION_OBJECT
    
    template<class T>
    struct constant
    {
        T val;
        constant(T v) : val(v) { }
        T operator()() const { return val; }
    };
    
    struct get_first
    {
        template<class T1, class T2>
        T1 operator()(std::pair<T1, T2> const & p) { return p.first; }
    };
    
    struct get_second
    {
        template<class T1, class T2>
        T2 operator()(std::pair<T1, T2> const & p) { return p.second; }
    };
    
} /* namespace */

#endif
