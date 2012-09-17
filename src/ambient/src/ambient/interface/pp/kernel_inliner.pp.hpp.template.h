#include <boost/preprocessor.hpp>
#define ARGS_MAX_LEN 11

#define cleanup_object(z, n, unused)                                                                                 \
    info<T ## n>::typed::deallocate<n>(o);

#define creditup_object(z, n, unused)                                                                                \
    info<T ## n>::typed::weight<n>(o); 

#define place_revision(z, n, unused)                                                                                 \
    info<T ## n>::typed::place<n>(o);                                                           

#define ready_revision(z, n, unused)                                                                                 \
    info<T ## n>::typed::ready<n>(o,e) &&

#define match_revision(z, n, unused)                                                                                 \
    info<T ## n>::typed::match<n>(o,t) ||

#define tag_revision(z, n, unused)                                                                                   \
    info<T ## n>::typed::tag<n>(o,t);

#define extract_arguments(z, n, unused)                                                                              \
    info<T ## n>::typed::modify<n>(arg ## n, o);

#define type_arg_list(z, n, pn)                                                                                      \
    BOOST_PP_COMMA_IF(n)                                                                                             \
    T ## n& arg ## n                                                                                                          

#define type_list(z, n, pn)                                                                                          \
    BOOST_PP_COMMA_IF(n)                                                                                             \
    T ## n&                                                                                                          

#define arg_list(z, n, pn)                                                                                           \
    BOOST_PP_COMMA_IF(n)                                                                                             \
    info<T ## n>::typed::revised<n>(o)
                                                                                                                     
#ifndef BOOST_PP_IS_ITERATING
#ifndef AMBIENT_INTERFACE_KERNEL_INLINER_PP
#define AMBIENT_INTERFACE_KERNEL_INLINER_PP


#define BOOST_PP_ITERATION_LIMITS (1, ARGS_MAX_LEN)
#define BOOST_PP_FILENAME_1 "ambient/interface/pp/kernel_inliner.pp.hpp.template.h"
#include BOOST_PP_ITERATE()
#endif
#else
#define n BOOST_PP_ITERATION()
#define TYPES_NUMBER n

template< class K, BOOST_PP_ENUM_PARAMS(TYPES_NUMBER, typename T) , void(K::*fp)( BOOST_PP_REPEAT(TYPES_NUMBER, type_list, BOOST_PP_ADD(n,1)) )>
struct kernel_inliner<void(K::*)( BOOST_PP_REPEAT(TYPES_NUMBER, type_list, BOOST_PP_ADD(n,1)) ), fp> {
    
    static inline void latch(cfunctor* o, BOOST_PP_REPEAT(TYPES_NUMBER, type_arg_list, n) ){
        BOOST_PP_REPEAT(TYPES_NUMBER, extract_arguments, ~) 
    }
    static inline void invoke(K* o){
        (o->*fp)( BOOST_PP_REPEAT(TYPES_NUMBER, arg_list, BOOST_PP_ADD(n,1)) );
    }
    static inline void cleanup(sfunctor* o){
        BOOST_PP_REPEAT(TYPES_NUMBER, cleanup_object, ~) 
    }
    static inline void weight(cfunctor* o){
        BOOST_PP_REPEAT(TYPES_NUMBER, creditup_object, ~) 
    }
    static inline void place(sfunctor* o){
        BOOST_PP_REPEAT(TYPES_NUMBER, place_revision, ~) 
    }
    static inline bool ready(sfunctor* o, void* e){
        return (BOOST_PP_REPEAT(TYPES_NUMBER, ready_revision, ~) true);
    }
    static inline bool match(sfunctor* o, void* t){
        return (BOOST_PP_REPEAT(TYPES_NUMBER, match_revision, ~) false);
    }
    static inline void tag(sfunctor* o, void* t){
        BOOST_PP_REPEAT(TYPES_NUMBER, tag_revision, ~)
    }
};

#endif
