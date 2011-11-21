#include <boost/preprocessor.hpp>
#define ARGS_MAX_LEN 10

#define type_list(z, n, pn)                                                                 \
    BOOST_PP_COMMA_IF(n)                                                                    \
    BOOST_PP_IF(BOOST_PP_EQUAL(n,pn), pinned,)                                              \
    T ## n&                                      

#define body_tn(z, n, text)                                                                 \
template < BOOST_PP_ENUM_PARAMS(TYPES_NUMBER, typename T) >                                 \
void extract_template(void (*)( BOOST_PP_REPEAT(TYPES_NUMBER, type_list, n) ));             \
template < BOOST_PP_ENUM_PARAMS(TYPES_NUMBER, typename T) >                                 \
void prototype_template(void (*)( BOOST_PP_REPEAT(TYPES_NUMBER, type_list, n) ));

#ifndef BOOST_PP_IS_ITERATING
#ifndef CONVERTOBJECTS_HPP
#define CONVERTOBJECTS_HPP
 
#define BOOST_PP_ITERATION_LIMITS (1, ARGS_MAX_LEN)
#define BOOST_PP_FILENAME_1 "ambient/core/operation/operation.pp.h.template.h"
#include BOOST_PP_ITERATE()
#endif
#else
#define n BOOST_PP_ITERATION()
#define TYPES_NUMBER n

template< typename FP, BOOST_PP_ENUM_PARAMS(TYPES_NUMBER, typename T) >
operation( FP op, BOOST_PP_ENUM_BINARY_PARAMS(TYPES_NUMBER, T, *arg) );
template < BOOST_PP_ENUM_PARAMS(TYPES_NUMBER, typename T) >     // specially for unpinned version
void extract_template(void (*)( BOOST_PP_REPEAT(TYPES_NUMBER, type_list, BOOST_PP_ADD(n,1)) ));
template < BOOST_PP_ENUM_PARAMS(TYPES_NUMBER, typename T) >     // specially for unpinned version
void prototype_template(void (*)( BOOST_PP_REPEAT(TYPES_NUMBER, type_list, BOOST_PP_ADD(n,1)) ));
BOOST_PP_REPEAT(TYPES_NUMBER, body_tn, ~) 

#endif
