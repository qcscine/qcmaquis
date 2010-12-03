#ifndef FUNCTION_OBJECTS_H
#define FUNCTION_OBJECTS_H

namespace functors {
    float conj(float v) { return v; }
    double conj(double v) { return v; }
    
#define DEFINE_FUNCTION_OBJECT(name, return_type, arg_type) \
struct f##name { template<class T> return_type operator() (arg_type t) { return name(t); } };
    
    DEFINE_FUNCTION_OBJECT(trace, typename T::value_type, T const &)
    DEFINE_FUNCTION_OBJECT(transpose, T, T const &)
    DEFINE_FUNCTION_OBJECT(conj, T, T)
    DEFINE_FUNCTION_OBJECT(conjugate, T, T const &)
    
#undef DEFINE_FUNCTION_OBJECT
    
    template<class T>
    struct constant
    {
        T val;
        constant(T v) : val(v) { }
        T operator()() const { return val; }
    };
    
} /* namespace */
    
#endif
