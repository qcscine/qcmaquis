#include <boost/preprocessor.hpp>
#define ARGS_MAX_LEN 10

#define arg_list(z, n, pn)                                                                                           \
    BOOST_PP_COMMA_IF(n)                                                                                             \
    info<T ## n>::unfold(arg ## n)

#ifndef BOOST_PP_IS_ITERATING
#ifndef AMBIENT_INTERFACE_PUSH_PP
#define AMBIENT_INTERFACE_PUSH_PP

namespace ambient{

template <typename K, class T0>
inline void push(T0& arg0){
    kernel_inliner<typename K::F,&K::c>::latch(new K(), info<T0>::unfold(arg0)); 
}
template <typename K, class T0, class T1>
inline void push(T0& arg0, T1& arg1){
    kernel_inliner<typename K::F,&K::c>::latch(new K(), info<T0>::unfold(arg0), info<T1>::unfold(arg1)); 
}
template <typename K, class T0, class T1, class T2>
inline void push(T0& arg0, T1& arg1, T2& arg2){
    kernel_inliner<typename K::F,&K::c>::latch(new K(), info<T0>::unfold(arg0), info<T1>::unfold(arg1), info<T2>::unfold(arg2)); 
}

#define BOOST_PP_ITERATION_LIMITS (4, ARGS_MAX_LEN)
#define BOOST_PP_FILENAME_1 "ambient/interface/pp/push.pp.hpp.template.h"
#include BOOST_PP_ITERATE()

} // namespace ambient
#endif
#else
#define n BOOST_PP_ITERATION()
#define TYPES_NUMBER n

template < typename K, BOOST_PP_ENUM_PARAMS(TYPES_NUMBER, class T) >
inline void push( BOOST_PP_ENUM_BINARY_PARAMS(TYPES_NUMBER, T, &arg) ){
    kernel_inliner<typename K::F,&K::c>::latch(new K(), BOOST_PP_REPEAT(TYPES_NUMBER, arg_list, BOOST_PP_ADD(n,1)) );
}

#endif
