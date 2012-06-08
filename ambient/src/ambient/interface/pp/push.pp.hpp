namespace ambient{

template <typename K, class T0>
inline void push(T0& arg0){
    kernel_inliner<typename K::F,K::c>::latch(new K(), info<T0>::unfold(arg0));
    ambient::playout();
}
template <typename K, class T0, class T1>
inline void push(T0& arg0, T1& arg1){
    kernel_inliner<typename K::F,K::c>::latch(new K(), info<T0>::unfold(arg0), info<T1>::unfold(arg1));
    ambient::playout();
}
template <typename K, class T0, class T1, class T2>
inline void push(T0& arg0, T1& arg1, T2& arg2){
    kernel_inliner<typename K::F,K::c>::latch(new K(), info<T0>::unfold(arg0), info<T1>::unfold(arg1), info<T2>::unfold(arg2));
    ambient::playout();
}
template < typename K, class T0 , class T1 , class T2 , class T3 >
inline void push(T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 ){
    kernel_inliner<typename K::F,K::c>::latch(new K(), info<T0>::unfold(arg0) , info<T1>::unfold(arg1) , info<T2>::unfold(arg2) , info<T3>::unfold(arg3) );
    ambient::playout();
}
template < typename K, class T0 , class T1 , class T2 , class T3 , class T4 >
inline void push(T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 , T4 &arg4 ){
    kernel_inliner<typename K::F,K::c>::latch(new K(), info<T0>::unfold(arg0) , info<T1>::unfold(arg1) , info<T2>::unfold(arg2) , info<T3>::unfold(arg3) , info<T4>::unfold(arg4) );
    ambient::playout();
}
template < typename K, class T0 , class T1 , class T2 , class T3 , class T4 , class T5 >
inline void push(T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 , T4 &arg4 , T5 &arg5 ){
    kernel_inliner<typename K::F,K::c>::latch(new K(), info<T0>::unfold(arg0) , info<T1>::unfold(arg1) , info<T2>::unfold(arg2) , info<T3>::unfold(arg3) , info<T4>::unfold(arg4) , info<T5>::unfold(arg5) );
    ambient::playout();
}
template < typename K, class T0 , class T1 , class T2 , class T3 , class T4 , class T5 , class T6 >
inline void push(T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 , T4 &arg4 , T5 &arg5 , T6 &arg6 ){
    kernel_inliner<typename K::F,K::c>::latch(new K(), info<T0>::unfold(arg0) , info<T1>::unfold(arg1) , info<T2>::unfold(arg2) , info<T3>::unfold(arg3) , info<T4>::unfold(arg4) , info<T5>::unfold(arg5) , info<T6>::unfold(arg6) );
    ambient::playout();
}

template < typename K, class T0 , class T1 , class T2 , class T3 , class T4 , class T5 , class T6 , class T7 >
inline void push(T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 , T4 &arg4 , T5 &arg5 , T6 &arg6 , T7 &arg7 ){
    kernel_inliner<typename K::F,K::c>::latch(new K(), info<T0>::unfold(arg0) , info<T1>::unfold(arg1) , info<T2>::unfold(arg2) , info<T3>::unfold(arg3) , info<T4>::unfold(arg4) , info<T5>::unfold(arg5) , info<T6>::unfold(arg6) , info<T7>::unfold(arg7) );
    ambient::playout();
}
template < typename K, class T0 , class T1 , class T2 , class T3 , class T4 , class T5 , class T6 , class T7 , class T8 >
inline void push(T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 , T4 &arg4 , T5 &arg5 , T6 &arg6 , T7 &arg7 , T8 &arg8 ){
    kernel_inliner<typename K::F,K::c>::latch(new K(), info<T0>::unfold(arg0) , info<T1>::unfold(arg1) , info<T2>::unfold(arg2) , info<T3>::unfold(arg3) , info<T4>::unfold(arg4) , info<T5>::unfold(arg5) , info<T6>::unfold(arg6) , info<T7>::unfold(arg7) , info<T8>::unfold(arg8) );
    ambient::playout();
}
template < typename K, class T0 , class T1 , class T2 , class T3 , class T4 , class T5 , class T6 , class T7 , class T8 , class T9 >
inline void push(T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 , T4 &arg4 , T5 &arg5 , T6 &arg6 , T7 &arg7 , T8 &arg8 , T9 &arg9 ){
    kernel_inliner<typename K::F,K::c>::latch(new K(), info<T0>::unfold(arg0) , info<T1>::unfold(arg1) , info<T2>::unfold(arg2) , info<T3>::unfold(arg3) , info<T4>::unfold(arg4) , info<T5>::unfold(arg5) , info<T6>::unfold(arg6) , info<T7>::unfold(arg7) , info<T8>::unfold(arg8) , info<T9>::unfold(arg9) );
    ambient::playout();
}

}
