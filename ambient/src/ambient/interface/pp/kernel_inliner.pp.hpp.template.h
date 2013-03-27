#include <boost/preprocessor.hpp>
#define ARGS_MAX_LEN 11

#define cleanup_object(z, n, unused)                                                                                 \
    info<T ## n>::typed::template deallocate<n>(o);

#define ready_revision(z, n, unused)                                                                                 \
    info<T ## n>::typed::template ready<n>(o) &&

#define extract_arguments(z, n, unused)                                                                              \
    info<T ## n>::typed::template modify<n>(arg ## n, o);

#define extract_local_arguments(z, n, unused)                                                                        \
    info<T ## n>::typed::template modify_local<n>(arg ## n, o);

#define extract_remote_arguments(z, n, unused)                                                                       \
    info<T ## n>::typed::template modify_remote<n>(arg ## n);

#define traverse_arguments(z, n, unused)                                                                             \
    info<T ## n>::typed::template pin<n>(o) ||

#define score_arguments(z, n, unused)                                                                                \
    info<T ## n>::typed::template score<n>(arg ## n);

#define type_arg_list(z, n, pn)                                                                                      \
    BOOST_PP_COMMA_IF(n)                                                                                             \
    T ## n& arg ## n                                                                                                          

#define type_list(z, n, pn)                                                                                          \
    BOOST_PP_COMMA_IF(n)                                                                                             \
    T ## n&                                                                                                          

#define arg_list(z, n, pn)                                                                                           \
    BOOST_PP_COMMA_IF(n)                                                                                             \
    info<T ## n>::typed::template revised<n>(o)
                                                                                                                     
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

template<BOOST_PP_ENUM_PARAMS(TYPES_NUMBER, typename T) , void(*fp)( BOOST_PP_REPEAT(TYPES_NUMBER, type_list, BOOST_PP_ADD(n,1)) )>
struct kernel_inliner<void(*)( BOOST_PP_REPEAT(TYPES_NUMBER, type_list, BOOST_PP_ADD(n,1)) ), fp> {
    static const int arity = n; 
    template<complexity O>
    static inline void latch(cfunctor* o, BOOST_PP_REPEAT(TYPES_NUMBER, type_arg_list, n) ){

        BOOST_PP_REPEAT(TYPES_NUMBER, score_arguments, ~)
        ambient::controller.schedule<O>();

        if(ambient::controller.remote()){
            BOOST_PP_REPEAT(TYPES_NUMBER, extract_remote_arguments, ~)
            return;
        }else if(ambient::controller.local()){
            BOOST_PP_REPEAT(TYPES_NUMBER, extract_local_arguments, ~) 
        }else{
            BOOST_PP_REPEAT(TYPES_NUMBER, extract_arguments, ~) 
        }

        BOOST_PP_REPEAT(TYPES_NUMBER, traverse_arguments, ~)
        ambient::controller.queue(o);
    }
    static inline void invoke(cfunctor* o){
        (*fp)( BOOST_PP_REPEAT(TYPES_NUMBER, arg_list, BOOST_PP_ADD(n,1)) );
    }
    static inline void cleanup(cfunctor* o){
        BOOST_PP_REPEAT(TYPES_NUMBER, cleanup_object, ~) 
    }
    static inline bool ready(cfunctor* o){
        return (BOOST_PP_REPEAT(TYPES_NUMBER, ready_revision, ~) true);
    }
};

#endif
