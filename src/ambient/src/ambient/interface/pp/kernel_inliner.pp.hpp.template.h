#include <boost/preprocessor.hpp>
#define ARGS_MAX_LEN 10

#define pin_object(z, n, pn)                                                                                         \
    BOOST_PP_IF(BOOST_PP_EQUAL(n,pn), o->pin = &ui_m_current(info<T ## n>::typed::dereference(o->arguments[n]));,)   \

#define cleanup_object(z, n, unused)                                                                                 \
    info<T ## n>::typed::deallocate(o->arguments[n]);

#define creditup_object(z, n, unused)                                                                                \
    info<T ## n>::typed::weight(o->arguments[n], o); 

#define place_revision(z, n, unused)                                                                                 \
    info<T ## n>::typed::place(o->arguments[n], o);                                                           

#define extract_arguments(z, n, unused)                                                                              \
    o->arguments[n] = (void*)info<T ## n>::typed::pointer(arg ## n);                                                 \
    o->revisions[n] = info<T ## n>::typed::modify(const_cast<T ## n&>(info<T ## n>::typed::dereference(o->arguments[n])), o);

#define type_arg_list(z, n, pn)                                                                                      \
    BOOST_PP_COMMA_IF(n)                                                                                             \
    T ## n& arg ## n                                                                                                          

#define type_list(z, n, pn)                                                                                          \
    BOOST_PP_COMMA_IF(n)                                                                                             \
    BOOST_PP_IF(BOOST_PP_EQUAL(n,pn), pinned,)                                                                       \
    T ## n&                                                                                                          

#define arg_list(z, n, pn)                                                                                           \
    BOOST_PP_COMMA_IF(n)                                                                                             \
    BOOST_PP_IF(BOOST_PP_EQUAL(n,pn), marked,)                                                                       \
    info<T ## n>::typed::revised(o->arguments[n], o->revisions[n])
                                                                                                                     
#define body_tn(z, n, text)                                                                                          \
template< BOOST_PP_ENUM_PARAMS(TYPES_NUMBER, typename T) , void(*fp)( BOOST_PP_REPEAT(TYPES_NUMBER, type_list, n) )> \
struct kernel_inliner<void(*)( BOOST_PP_REPEAT(TYPES_NUMBER, type_list, n) ), fp> {                               \
                                                                                                                     \
    static inline void latch(cfunctor* o, BOOST_PP_REPEAT(TYPES_NUMBER, type_arg_list, n) ){                         \
        o->arguments    = (void**)malloc(sizeof(void*)*TYPES_NUMBER);                                                \
        o->revisions    = (size_t*)malloc(sizeof(size_t)*TYPES_NUMBER);                                              \
        BOOST_PP_REPEAT(TYPES_NUMBER, extract_arguments, ~)                                                          \
        BOOST_PP_REPEAT(TYPES_NUMBER, pin_object, n)                                                                 \
        ambient::controller.push(o);                                                                                 \
    }                                                                                                                \
    static inline void invoke(sfunctor* o){                                                                          \
        fp( BOOST_PP_REPEAT(TYPES_NUMBER, arg_list, n) );                                                            \
    }                                                                                                                \
    static inline void cleanup(sfunctor* o){                                                                         \
        BOOST_PP_REPEAT(TYPES_NUMBER, cleanup_object, ~)                                                             \
    }                                                                                                                \
    static inline void weight(cfunctor* o){                                                                          \
        BOOST_PP_REPEAT(TYPES_NUMBER, creditup_object, ~)                                                            \
    }                                                                                                                \
    static inline void place(sfunctor* o){                                                                           \
        BOOST_PP_REPEAT(TYPES_NUMBER, place_revision, ~)                                                             \
    }                                                                                                                \
    static inline bool pretend(cfunctor* o){                                                                         \
        return false;                                                                                                \
    }                                                                                                                \
    static inline void l_check_complete(cfunctor* o){                                                                \
        o->check_complete();                                                                                         \
    }                                                                                                                \
    static inline void c_check_complete(cfunctor* o){                                                                \
        o->check_complete();                                                                                         \
    }                                                                                                                \
                                                                                                                     \
};


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

template< BOOST_PP_ENUM_PARAMS(TYPES_NUMBER, typename T) , void(*fp)( BOOST_PP_REPEAT(TYPES_NUMBER, type_list, BOOST_PP_ADD(n,1)) )>
struct kernel_inliner<void(*)( BOOST_PP_REPEAT(TYPES_NUMBER, type_list, BOOST_PP_ADD(n,1)) ), fp> {

    static inline void latch(cfunctor* o, BOOST_PP_REPEAT(TYPES_NUMBER, type_arg_list, n) ){
        o->arguments    = (void**)malloc(sizeof(void*)*TYPES_NUMBER);
        o->revisions    = (size_t*)malloc(sizeof(size_t)*TYPES_NUMBER);
        BOOST_PP_REPEAT(TYPES_NUMBER, extract_arguments, ~) 
        o->pin = NULL;
        ambient::controller.push(o); 
    }
    static inline void invoke(sfunctor* o){
        fp( BOOST_PP_REPEAT(TYPES_NUMBER, arg_list, BOOST_PP_ADD(n,1)) );
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
    static inline bool pretend(cfunctor* o){
        o->lock();
        o->workload--;
        o->unlock();
        if(!o->workload) return false;
        return true;
    }
    static inline void l_check_complete(cfunctor* o){
        o->lock();
        if(o->workload == 1) controller.execute_free_mod(o);
        else o->workload--;
        o->unlock();
    }
    static inline void c_check_complete(cfunctor* o){
        o->check_complete();
    }

};

BOOST_PP_REPEAT(n, body_tn, ~) 

#endif
