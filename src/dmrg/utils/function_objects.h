/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef FUNCTION_OBJECTS_H
#define FUNCTION_OBJECTS_H

namespace utils {
    inline float conj(float v) { return v; }
    inline double conj(double v) { return v; }

#define DEFINE_FUNCTION_OBJECT(name, return_type, arg_type) \
struct functor_##name { template<class T> return_type operator() (arg_type t) { return name(t); } };
    
    DEFINE_FUNCTION_OBJECT(trace, typename T::value_type, T const &)
    DEFINE_FUNCTION_OBJECT(transpose, T, T const &)
    DEFINE_FUNCTION_OBJECT(conj, T, T)
    DEFINE_FUNCTION_OBJECT(conjugate, T, T const &)
    DEFINE_FUNCTION_OBJECT(sqrt, T, T const &)
    DEFINE_FUNCTION_OBJECT(exp, T, T const &)
    DEFINE_FUNCTION_OBJECT(size_of, std::size_t, T const &)
    
#undef DEFINE_FUNCTION_OBJECT
    
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
