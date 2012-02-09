#include <boost/preprocessor.hpp>
#define ARGS_MAX_LEN 10

#ifndef BOOST_PP_IS_ITERATING
#ifndef CONVERTOBJECTS_HPP
#define CONVERTOBJECTS_HPP

template <typename FP, class T0>
void push(FP l_kernel, FP c_kernel, T0& arg0){
    ambient::controller.push(new models::operation(l_kernel, c_kernel, arg0)); 
}
template <typename FP, class T0, class T1>
void push(FP l_kernel, FP c_kernel, T0& arg0, T1& arg1){
    ambient::controller.push(new models::operation(l_kernel, c_kernel, arg0, arg1));
}
template <typename FP, class T0, class T1, class T2>
void push(FP l_kernel, FP c_kernel, T0& arg0, T1& arg1, T2& arg2){
    ambient::controller.push(new models::operation(l_kernel, c_kernel, arg0, arg1, arg2)); 
}

#define BOOST_PP_ITERATION_LIMITS (4, ARGS_MAX_LEN)
#define BOOST_PP_FILENAME_1 "ambient/interface/pp/push.pp.hpp.template.h"
#include BOOST_PP_ITERATE()
#endif
#else
#define n BOOST_PP_ITERATION()
#define TYPES_NUMBER n

template < typename FP, BOOST_PP_ENUM_PARAMS(TYPES_NUMBER, class T) >
void push( FP l_kernel, FP c_kernel, BOOST_PP_ENUM_BINARY_PARAMS(TYPES_NUMBER, T, &arg) ){
    ambient::controller.push(new models::operation( l_kernel, c_kernel, BOOST_PP_ENUM_PARAMS(TYPES_NUMBER, arg) ));
}

#endif
