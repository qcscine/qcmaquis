#include <boost/preprocessor.hpp>
#define ARGS_MAX_LEN 4

#define extract_arguments(z, n, unused)                                                     \
    this->arguments[n] = boost::shared_ptr<void>(arg ## n);

#define type_list(z, n, pn)                                                                 \
    BOOST_PP_COMMA_IF(n)                                                                    \
    BOOST_PP_IF(BOOST_PP_EQUAL(n,pn), pinned,)                                              \
    T ## n&                                      

#define arg_list(z, n, pn)                                                                  \
    BOOST_PP_COMMA_IF(n)                                                                    \
    BOOST_PP_IF(BOOST_PP_EQUAL(n,pn), marked,)                                              \
    *static_cast<T ## n*>(this->arguments[n].get())

#define body_tn(z, n, text)                                                                 \
template < BOOST_PP_ENUM_PARAMS(TYPES_NUMBER, typename T) >                                 \
void operation::prototype_template(void (*)( BOOST_PP_REPEAT(TYPES_NUMBER, type_list, n) )) \
{                                                                                           \
    ((void (*)( BOOST_PP_REPEAT(TYPES_NUMBER, type_list, n) ))this->operation_ptr)          \
    ( BOOST_PP_REPEAT(TYPES_NUMBER, arg_list, n) );                                         \
}

#ifndef BOOST_PP_IS_ITERATING
#ifndef CONVERTOBJECTS_HPP
#define CONVERTOBJECTS_HPP

namespace ambient{ namespace core{

#define BOOST_PP_ITERATION_LIMITS (1, ARGS_MAX_LEN)
#define BOOST_PP_FILENAME_1 "ambient/core/operation/operation.pp.sa.hpp.template.h"
#include BOOST_PP_ITERATE()
#endif
#else
#define n BOOST_PP_ITERATION()
#define TYPES_NUMBER n

template< typename FP, BOOST_PP_ENUM_PARAMS(TYPES_NUMBER, typename T) >
operation::operation( FP op, BOOST_PP_ENUM_BINARY_PARAMS(TYPES_NUMBER, T, *arg) ){
    this->init();
    this->operation_ptr = (void(*)())op;
    this->count = TYPES_NUMBER;
    this->profiles  = (p_object**)malloc(sizeof(p_object*)*this->count);
    this->arguments = new boost::shared_ptr<void>[this->count];
    BOOST_PP_REPEAT(TYPES_NUMBER, extract_arguments, ~) 
    void(operation::*ptr)(FP); ptr = &operation::prototype_template;
    this->prototype = (void(operation::*)())ptr;
}

template < BOOST_PP_ENUM_PARAMS(TYPES_NUMBER, typename T) >
void operation::prototype_template(void (*)( BOOST_PP_REPEAT(TYPES_NUMBER, type_list, BOOST_PP_ADD(n,1)) ))
{
    ((void (*)( BOOST_PP_REPEAT(TYPES_NUMBER, type_list, BOOST_PP_ADD(n,1)) ))this->operation_ptr)
    ( BOOST_PP_REPEAT(TYPES_NUMBER, arg_list, BOOST_PP_ADD(n,1)) );
}
BOOST_PP_REPEAT(n, body_tn, ~) 

BOOST_PP_IF(BOOST_PP_EQUAL(TYPES_NUMBER, ARGS_MAX_LEN), } },)
#endif
