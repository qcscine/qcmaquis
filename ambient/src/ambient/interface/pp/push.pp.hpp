template <typename FP, class T0>
void push(FP l_kernel, FP c_kernel, T0& arg0){
    ambient::controller.push(new models::operation(l_kernel, c_kernel, models::info<T0>::unfold(arg0)));
}
template <typename FP, class T0, class T1>
void push(FP l_kernel, FP c_kernel, T0& arg0, T1& arg1){
    ambient::controller.push(new models::operation(l_kernel, c_kernel, models::info<T0>::unfold(arg0), models::info<T1>::unfold(arg1)));
}
template <typename FP, class T0, class T1, class T2>
void push(FP l_kernel, FP c_kernel, T0& arg0, T1& arg1, T2& arg2){
    ambient::controller.push(new models::operation(l_kernel, c_kernel, models::info<T0>::unfold(arg0), models::info<T1>::unfold(arg1), models::info<T2>::unfold(arg2)));
}
template < typename FP, class T0 , class T1 , class T2 , class T3 >
void push( FP l_kernel, FP c_kernel, T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 ){
    ambient::controller.push(new models::operation( l_kernel, c_kernel, models::info<T0>::unfold(arg0) , models::info<T1>::unfold(arg1) , models::info<T2>::unfold(arg2) , models::info<T3>::unfold(arg3) ));
}
template < typename FP, class T0 , class T1 , class T2 , class T3 , class T4 >
void push( FP l_kernel, FP c_kernel, T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 , T4 &arg4 ){
    ambient::controller.push(new models::operation( l_kernel, c_kernel, models::info<T0>::unfold(arg0) , models::info<T1>::unfold(arg1) , models::info<T2>::unfold(arg2) , models::info<T3>::unfold(arg3) , models::info<T4>::unfold(arg4) ));
}
template < typename FP, class T0 , class T1 , class T2 , class T3 , class T4 , class T5 >
void push( FP l_kernel, FP c_kernel, T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 , T4 &arg4 , T5 &arg5 ){
    ambient::controller.push(new models::operation( l_kernel, c_kernel, models::info<T0>::unfold(arg0) , models::info<T1>::unfold(arg1) , models::info<T2>::unfold(arg2) , models::info<T3>::unfold(arg3) , models::info<T4>::unfold(arg4) , models::info<T5>::unfold(arg5) ));
}
template < typename FP, class T0 , class T1 , class T2 , class T3 , class T4 , class T5 , class T6 >
void push( FP l_kernel, FP c_kernel, T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 , T4 &arg4 , T5 &arg5 , T6 &arg6 ){
    ambient::controller.push(new models::operation( l_kernel, c_kernel, models::info<T0>::unfold(arg0) , models::info<T1>::unfold(arg1) , models::info<T2>::unfold(arg2) , models::info<T3>::unfold(arg3) , models::info<T4>::unfold(arg4) , models::info<T5>::unfold(arg5) , models::info<T6>::unfold(arg6) ));
}

template < typename FP, class T0 , class T1 , class T2 , class T3 , class T4 , class T5 , class T6 , class T7 >
void push( FP l_kernel, FP c_kernel, T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 , T4 &arg4 , T5 &arg5 , T6 &arg6 , T7 &arg7 ){
    ambient::controller.push(new models::operation( l_kernel, c_kernel, models::info<T0>::unfold(arg0) , models::info<T1>::unfold(arg1) , models::info<T2>::unfold(arg2) , models::info<T3>::unfold(arg3) , models::info<T4>::unfold(arg4) , models::info<T5>::unfold(arg5) , models::info<T6>::unfold(arg6) , models::info<T7>::unfold(arg7) ));
}
template < typename FP, class T0 , class T1 , class T2 , class T3 , class T4 , class T5 , class T6 , class T7 , class T8 >
void push( FP l_kernel, FP c_kernel, T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 , T4 &arg4 , T5 &arg5 , T6 &arg6 , T7 &arg7 , T8 &arg8 ){
    ambient::controller.push(new models::operation( l_kernel, c_kernel, models::info<T0>::unfold(arg0) , models::info<T1>::unfold(arg1) , models::info<T2>::unfold(arg2) , models::info<T3>::unfold(arg3) , models::info<T4>::unfold(arg4) , models::info<T5>::unfold(arg5) , models::info<T6>::unfold(arg6) , models::info<T7>::unfold(arg7) , models::info<T8>::unfold(arg8) ));
}
template < typename FP, class T0 , class T1 , class T2 , class T3 , class T4 , class T5 , class T6 , class T7 , class T8 , class T9 >
void push( FP l_kernel, FP c_kernel, T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 , T4 &arg4 , T5 &arg5 , T6 &arg6 , T7 &arg7 , T8 &arg8 , T9 &arg9 ){
    ambient::controller.push(new models::operation( l_kernel, c_kernel, models::info<T0>::unfold(arg0) , models::info<T1>::unfold(arg1) , models::info<T2>::unfold(arg2) , models::info<T3>::unfold(arg3) , models::info<T4>::unfold(arg4) , models::info<T5>::unfold(arg5) , models::info<T6>::unfold(arg6) , models::info<T7>::unfold(arg7) , models::info<T8>::unfold(arg8) , models::info<T9>::unfold(arg9) ));
}
