template <typename FP, class T0>
void push(FP l_kernel, FP c_kernel, T0& arg0){
    static __a_timer time("ambient_push_1"); time.begin();
    ambient::controller.push(new models::operation(l_kernel, c_kernel, models::info<T0>::unfold(arg0)));
    time.end();
    ambient::playout();
}
template <typename FP, class T0, class T1>
void push(FP l_kernel, FP c_kernel, T0& arg0, T1& arg1){
    static __a_timer time("ambient_push_2"); time.begin();
    ambient::controller.push(new models::operation(l_kernel, c_kernel, models::info<T0>::unfold(arg0), models::info<T1>::unfold(arg1)));
    time.end();
    ambient::playout();
}
template <typename FP, class T0, class T1, class T2>
void push(FP l_kernel, FP c_kernel, T0& arg0, T1& arg1, T2& arg2){
    static __a_timer time("ambient_push_3"); time.begin();
    ambient::controller.push(new models::operation(l_kernel, c_kernel, models::info<T0>::unfold(arg0), models::info<T1>::unfold(arg1), models::info<T2>::unfold(arg2)));
    time.end();
    ambient::playout();
}
template < typename FP, class T0 , class T1 , class T2 , class T3 >
void push( FP l_kernel, FP c_kernel, T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 ){
    static __a_timer time("ambient_push_4"); time.begin();
    ambient::controller.push(new models::operation( l_kernel, c_kernel, models::info<T0>::unfold(arg0) , models::info<T1>::unfold(arg1) , models::info<T2>::unfold(arg2) , models::info<T3>::unfold(arg3) ));
    time.end();
    ambient::playout();
}
template < typename FP, class T0 , class T1 , class T2 , class T3 , class T4 >
void push( FP l_kernel, FP c_kernel, T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 , T4 &arg4 ){
    static __a_timer time("ambient_push_5"); time.begin();
    ambient::controller.push(new models::operation( l_kernel, c_kernel, models::info<T0>::unfold(arg0) , models::info<T1>::unfold(arg1) , models::info<T2>::unfold(arg2) , models::info<T3>::unfold(arg3) , models::info<T4>::unfold(arg4) ));
    time.end();
    ambient::playout();
}
template < typename FP, class T0 , class T1 , class T2 , class T3 , class T4 , class T5 >
void push( FP l_kernel, FP c_kernel, T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 , T4 &arg4 , T5 &arg5 ){
    static __a_timer time("ambient_push_6"); time.begin();
    ambient::controller.push(new models::operation( l_kernel, c_kernel, models::info<T0>::unfold(arg0) , models::info<T1>::unfold(arg1) , models::info<T2>::unfold(arg2) , models::info<T3>::unfold(arg3) , models::info<T4>::unfold(arg4) , models::info<T5>::unfold(arg5) ));
    time.end();
    ambient::playout();
}
template < typename FP, class T0 , class T1 , class T2 , class T3 , class T4 , class T5 , class T6 >
void push( FP l_kernel, FP c_kernel, T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 , T4 &arg4 , T5 &arg5 , T6 &arg6 ){
    static __a_timer time("ambient_push_7"); time.begin();
    ambient::controller.push(new models::operation( l_kernel, c_kernel, models::info<T0>::unfold(arg0) , models::info<T1>::unfold(arg1) , models::info<T2>::unfold(arg2) , models::info<T3>::unfold(arg3) , models::info<T4>::unfold(arg4) , models::info<T5>::unfold(arg5) , models::info<T6>::unfold(arg6) ));
    time.end();
    ambient::playout();
}

template < typename FP, class T0 , class T1 , class T2 , class T3 , class T4 , class T5 , class T6 , class T7 >
void push( FP l_kernel, FP c_kernel, T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 , T4 &arg4 , T5 &arg5 , T6 &arg6 , T7 &arg7 ){
    static __a_timer time("ambient_push_8"); time.begin();
    ambient::controller.push(new models::operation( l_kernel, c_kernel, models::info<T0>::unfold(arg0) , models::info<T1>::unfold(arg1) , models::info<T2>::unfold(arg2) , models::info<T3>::unfold(arg3) , models::info<T4>::unfold(arg4) , models::info<T5>::unfold(arg5) , models::info<T6>::unfold(arg6) , models::info<T7>::unfold(arg7) ));
    time.end();
    ambient::playout();
}
template < typename FP, class T0 , class T1 , class T2 , class T3 , class T4 , class T5 , class T6 , class T7 , class T8 >
void push( FP l_kernel, FP c_kernel, T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 , T4 &arg4 , T5 &arg5 , T6 &arg6 , T7 &arg7 , T8 &arg8 ){
    static __a_timer time("ambient_push_9"); time.begin();
    ambient::controller.push(new models::operation( l_kernel, c_kernel, models::info<T0>::unfold(arg0) , models::info<T1>::unfold(arg1) , models::info<T2>::unfold(arg2) , models::info<T3>::unfold(arg3) , models::info<T4>::unfold(arg4) , models::info<T5>::unfold(arg5) , models::info<T6>::unfold(arg6) , models::info<T7>::unfold(arg7) , models::info<T8>::unfold(arg8) ));
    time.end();
    ambient::playout();
}
template < typename FP, class T0 , class T1 , class T2 , class T3 , class T4 , class T5 , class T6 , class T7 , class T8 , class T9 >
void push( FP l_kernel, FP c_kernel, T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 , T4 &arg4 , T5 &arg5 , T6 &arg6 , T7 &arg7 , T8 &arg8 , T9 &arg9 ){
    static __a_timer time("ambient_push_10"); time.begin();
    ambient::controller.push(new models::operation( l_kernel, c_kernel, models::info<T0>::unfold(arg0) , models::info<T1>::unfold(arg1) , models::info<T2>::unfold(arg2) , models::info<T3>::unfold(arg3) , models::info<T4>::unfold(arg4) , models::info<T5>::unfold(arg5) , models::info<T6>::unfold(arg6) , models::info<T7>::unfold(arg7) , models::info<T8>::unfold(arg8) , models::info<T9>::unfold(arg9) ));
    time.end();
    ambient::playout();
}
