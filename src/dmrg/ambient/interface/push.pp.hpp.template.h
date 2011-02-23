#include <boost/preprocessor.hpp>
#define ARGS_MAX_LEN 4

#define type_list(z, n, pn)                                                                 \
    BOOST_PP_COMMA_IF(n)                                                                    \
    BOOST_PP_IF(BOOST_PP_EQUAL(n,pn), pinned,)                                              \
    T ## n&                                      

#define body_tn(z, n, text)                                                                 \
template < BOOST_PP_ENUM_PARAMS(TYPES_NUMBER, typename T) >                                 \
void prototype_template(void (*)( BOOST_PP_REPEAT(TYPES_NUMBER, type_list, n) ));

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
template <typename FL, typename FC, class T0, class T1, class T2>
void push(FL l_kernel, FC c_kernel, T0& arg0, T1& arg1, T2& arg2){
    if(get_profile(arg1)->proxy) pin(arg0, arg2);
    if(get_profile(arg2)->proxy) pin(arg1, arg2);
    ambient::engine.push(new core::operation(l_kernel, &arg0, &arg1, &arg2), 
                         new core::operation(c_kernel, &arg0, &arg1, &arg2)); 
}
template <typename ST, typename FL, typename FC, class T0, class T1> 
ST push(FL l_kernel, FC c_kernel, T0& arg0, T1& arg1){
    void_pt* handle = new void_pt((ST*)NULL);
    ST out(handle);
    push(l_kernel, c_kernel, arg0, arg1, out);
    return out;
}

#define BOOST_PP_ITERATION_LIMITS (4, ARGS_MAX_LEN)
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
