#include <boost/preprocessor.hpp>
#define ARGS_MAX_LEN 10

#define pin_object(z, n, pn)                                                                \
    BOOST_PP_IF(BOOST_PP_EQUAL(n,pn), this->pin = this->arguments[n];,)

#define cleanup_object(z, n, unused)                                                        \
    info<T ## n>::typed::deallocate(this->arguments[n]);

#define creditup_object(z, n, unused)                                                       \
    info<T ## n>::typed::weight(this->arguments[n], this);

#define extract_arguments(z, n, unused)                                                     \
    this->arguments[n] = (void*)info<T ## n>::typed::pointer(arg ## n);                     \
    info<T ## n>::typed::modify(arg ## n, this);                                            \

#define type_list(z, n, pn)                                                                 \
    BOOST_PP_COMMA_IF(n)                                                                    \
    BOOST_PP_IF(BOOST_PP_EQUAL(n,pn), pinned,)                                              \
    T ## n&                                     

#define arg_list(z, n, pn)                                                                  \
    BOOST_PP_COMMA_IF(n)                                                                    \
    BOOST_PP_IF(BOOST_PP_EQUAL(n,pn), marked,)                                              \
    info<T ## n>::typed::dereference(this->arguments[n])

#define body_tn(z, n, text)                                                                 \
template < BOOST_PP_ENUM_PARAMS(TYPES_NUMBER, typename T) >                                 \
void prototype_template(void (*)( BOOST_PP_REPEAT(TYPES_NUMBER, type_list, n) )) \
{                                                                                           \
    ( (void (*)( BOOST_PP_REPEAT(TYPES_NUMBER, type_list, n) )) this->op )                  \
    ( BOOST_PP_REPEAT(TYPES_NUMBER, arg_list, n) );                                         \
}                                                                                           \
template < BOOST_PP_ENUM_PARAMS(TYPES_NUMBER, typename T) >                                 \
void cleanup_template(void (*)( BOOST_PP_REPEAT(TYPES_NUMBER, type_list, n) ))   \
{                                                                                           \
    BOOST_PP_REPEAT(TYPES_NUMBER, cleanup_object, ~)                                        \
}                                                                                           \
template < BOOST_PP_ENUM_PARAMS(TYPES_NUMBER, typename T) >                                 \
void creditup_template(void (*)( BOOST_PP_REPEAT(TYPES_NUMBER, type_list, n) ))   \
{                                                                                           \
    BOOST_PP_REPEAT(TYPES_NUMBER, creditup_object, ~)                                        \
}                                                                                           \
template < BOOST_PP_ENUM_PARAMS(TYPES_NUMBER, typename T) >                                 \
void mark_pin(void (*)( BOOST_PP_REPEAT(TYPES_NUMBER, type_list, n) ))           \
{                                                                                           \
    BOOST_PP_REPEAT(TYPES_NUMBER, pin_object, n)                                            \
}

#ifndef BOOST_PP_IS_ITERATING
#ifndef CONVERTOBJECTS_HPP
#define CONVERTOBJECTS_HPP

#define DEFINE_MARKED #define marked NULL,
#define DEFINE_PINNED #define pinned ambient::models::ambient_pin* ,
#define UNDEF_MARKED #undef marked
#define UNDEF_PINNED #undef pinned
DEFINE_MARKED
DEFINE_PINNED

#define BOOST_PP_ITERATION_LIMITS (1, ARGS_MAX_LEN)
#define BOOST_PP_FILENAME_1 "ambient/models/operation/pp/operation.pp.hpp.template.h"
#include BOOST_PP_ITERATE()
#endif
#else
#define n BOOST_PP_ITERATION()
#define TYPES_NUMBER n

template< typename FP, BOOST_PP_ENUM_PARAMS(TYPES_NUMBER, typename T) >
operation( FP logistics, FP computing, BOOST_PP_ENUM_BINARY_PARAMS(TYPES_NUMBER, T, &arg) ){
    this->logistics_ptr = (void(*)())logistics;
    this->computing_ptr = (void(*)())computing;
    this->op            = this->logistics_ptr;
    this->credit    = 0;
    this->state     = MARKUP;
    this->count     = TYPES_NUMBER;
    this->arguments = (void**)malloc(sizeof(void*)*this->count);
    BOOST_PP_REPEAT(TYPES_NUMBER, extract_arguments, ~) 
    void(operation::*ptr)(FP); ptr = &operation::prototype_template;
    this->prototype = (void(operation::*)())ptr;
    ptr = &operation::cleanup_template;
    this->cleanup = (void(operation::*)())ptr;
    ptr = &operation::creditup_template;
    this->creditup = (void(operation::*)())ptr;
    this->mark_pin(logistics);
}

template < BOOST_PP_ENUM_PARAMS(TYPES_NUMBER, typename T) >
void prototype_template(void (*)( BOOST_PP_REPEAT(TYPES_NUMBER, type_list, BOOST_PP_ADD(n,1)) )){
    ( (void (*)( BOOST_PP_REPEAT(TYPES_NUMBER, type_list, BOOST_PP_ADD(n,1)) )) this->op )
    ( BOOST_PP_REPEAT(TYPES_NUMBER, arg_list, BOOST_PP_ADD(n,1)) );
}

template < BOOST_PP_ENUM_PARAMS(TYPES_NUMBER, typename T) >
void cleanup_template(void (*)( BOOST_PP_REPEAT(TYPES_NUMBER, type_list, BOOST_PP_ADD(n,1)) )){
    BOOST_PP_REPEAT(TYPES_NUMBER, cleanup_object, ~) 
}

template < BOOST_PP_ENUM_PARAMS(TYPES_NUMBER, typename T) >
void creditup_template(void (*)( BOOST_PP_REPEAT(TYPES_NUMBER, type_list, BOOST_PP_ADD(n,1)) )){
    BOOST_PP_REPEAT(TYPES_NUMBER, creditup_object, ~) 
}

template < BOOST_PP_ENUM_PARAMS(TYPES_NUMBER, typename T) >
void mark_pin(void (*)( BOOST_PP_REPEAT(TYPES_NUMBER, type_list, BOOST_PP_ADD(n,1)) )){
}
BOOST_PP_REPEAT(n, body_tn, ~) 

#endif

#ifndef BOOST_PP_IS_ITERATING
UNDEF_MARKED
UNDEF_PINNED
#endif
