template <typename FL, typename FC, class T0>
void push(FL l_kernel, FC c_kernel, T0& arg0){
    ambient::engine.push(new core::operation(l_kernel, &arg0),
                         new core::operation(c_kernel, &arg0));
}
template <typename FL, typename FC, class T0, class T1>
void push(FL l_kernel, FC c_kernel, T0& arg0, T1& arg1){
    ambient::engine.push(new core::operation(l_kernel, &arg0, &arg1),
                         new core::operation(c_kernel, &arg0, &arg1));
}
template <typename FL, typename FC, class T0, class T1, class T2>
void push(FL l_kernel, FC c_kernel, T0& arg0, T1& arg1, T2& arg2){
    ambient::engine.push(new core::operation(l_kernel, &arg0, &arg1, &arg2),
			 new core::operation(c_kernel, &arg0, &arg1, &arg2));
}
template < typename FL, typename FC, class T0 , class T1 , class T2 , class T3 >
void push( FL l_kernel, FC c_kernel, T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 ){
    ambient::engine.push(new core::operation( l_kernel, &arg0 , &arg1 , &arg2 , &arg3 ),
                         new core::operation( c_kernel, &arg0 , &arg1 , &arg2 , &arg3 ));
}
template < typename FL, typename FC, class T0 , class T1 , class T2 , class T3 , class T4 >
void push( FL l_kernel, FC c_kernel, T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 , T4 &arg4 ){
    ambient::engine.push(new core::operation( l_kernel, &arg0 , &arg1 , &arg2 , &arg3 , &arg4 ),
                         new core::operation( c_kernel, &arg0 , &arg1 , &arg2 , &arg3 , &arg4 ));
}
template < typename FL, typename FC, class T0 , class T1 , class T2 , class T3 , class T4 , class T5 >
void push( FL l_kernel, FC c_kernel, T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 , T4 &arg4 , T5 &arg5 ){
    ambient::engine.push(new core::operation( l_kernel, &arg0 , &arg1 , &arg2 , &arg3 , &arg4 , &arg5 ),
                         new core::operation( c_kernel, &arg0 , &arg1 , &arg2 , &arg3 , &arg4 , &arg5 ));
}



template < typename FL, typename FC, class T0 , class T1 , class T2 , class T3 , class T4 , class T5 , class T6 >
void push( FL l_kernel, FC c_kernel, T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 , T4 &arg4 , T5 &arg5 , T6 &arg6 ){
    ambient::engine.push(new core::operation( l_kernel, &arg0 , &arg1 , &arg2 , &arg3 , &arg4 , &arg5 , &arg6 ),
                         new core::operation( c_kernel, &arg0 , &arg1 , &arg2 , &arg3 , &arg4 , &arg5 , &arg6 ));
}
template < typename FL, typename FC, class T0 , class T1 , class T2 , class T3 , class T4 , class T5 , class T6 , class T7 >
void push( FL l_kernel, FC c_kernel, T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 , T4 &arg4 , T5 &arg5 , T6 &arg6 , T7 &arg7 ){
    ambient::engine.push(new core::operation( l_kernel, &arg0 , &arg1 , &arg2 , &arg3 , &arg4 , &arg5 , &arg6 , &arg7 ),
                         new core::operation( c_kernel, &arg0 , &arg1 , &arg2 , &arg3 , &arg4 , &arg5 , &arg6 , &arg7 ));
}
template < typename FL, typename FC, class T0 , class T1 , class T2 , class T3 , class T4 , class T5 , class T6 , class T7 , class T8 >
void push( FL l_kernel, FC c_kernel, T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 , T4 &arg4 , T5 &arg5 , T6 &arg6 , T7 &arg7 , T8 &arg8 ){
    ambient::engine.push(new core::operation( l_kernel, &arg0 , &arg1 , &arg2 , &arg3 , &arg4 , &arg5 , &arg6 , &arg7 , &arg8 ),
                         new core::operation( c_kernel, &arg0 , &arg1 , &arg2 , &arg3 , &arg4 , &arg5 , &arg6 , &arg7 , &arg8 ));
}
template < typename FL, typename FC, class T0 , class T1 , class T2 , class T3 , class T4 , class T5 , class T6 , class T7 , class T8 , class T9 >
void push( FL l_kernel, FC c_kernel, T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 , T4 &arg4 , T5 &arg5 , T6 &arg6 , T7 &arg7 , T8 &arg8 , T9 &arg9 ){
    ambient::engine.push(new core::operation( l_kernel, &arg0 , &arg1 , &arg2 , &arg3 , &arg4 , &arg5 , &arg6 , &arg7 , &arg8 , &arg9 ),
                         new core::operation( c_kernel, &arg0 , &arg1 , &arg2 , &arg3 , &arg4 , &arg5 , &arg6 , &arg7 , &arg8 , &arg9 ));
}
