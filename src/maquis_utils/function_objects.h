/*****************************************************************************
 *
 * MAQUIS MATRICES Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
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
