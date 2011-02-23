#include <boost/preprocessor.hpp>
#define ARGS_MAX_LEN 4

#ifndef BOOST_PP_IS_ITERATING
#ifndef CONVERTOBJECTS_HPP
#define CONVERTOBJECTS_HPP

template <typename FL, typename FC, class T0>
void push(FL l_kernel, FC c_kernel, T0& arg0){
    ambient::engine.push(new core::operation(l_kernel, &arg0), 
                         new core::operation(c_kernel, &arg0)); 
}
template <typename FL, typename FC, class T0, class T1>
void push(FL l_kernel, FC c_kernel, T0& arg0, T1& arg1){
    ambient::engine.push(new core::operation(l_kernel, &arg0, &arg1), 
                         new core::operation(c_kernel, &arg0, &arg1)); 
}

#define BOOST_PP_ITERATION_LIMITS (4, ARGS_MAX_LEN) // 3 is special case
#define BOOST_PP_FILENAME_1 "ambient/interface/push.pp.hpp.template.h"
#include BOOST_PP_ITERATE()
#endif
#else
#define n BOOST_PP_ITERATION()
#define TYPES_NUMBER n

template < typename FL, typename FC, BOOST_PP_ENUM_PARAMS(TYPES_NUMBER, class T) >
void push( FL l_kernel, FC c_kernel, BOOST_PP_ENUM_BINARY_PARAMS(TYPES_NUMBER, T, &arg) ){
    ambient::engine.push(new core::operation( l_kernel, BOOST_PP_ENUM_PARAMS(TYPES_NUMBER, &arg) ), 
                         new core::operation( c_kernel, BOOST_PP_ENUM_PARAMS(TYPES_NUMBER, &arg) )); 
}

#endif
