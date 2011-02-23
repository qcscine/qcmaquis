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
template < typename FL, typename FC, class T0 , class T1 , class T2 , class T3 >
void push( FL l_kernel, FC c_kernel, T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 ){
    ambient::engine.push(new core::operation( l_kernel, &arg0 , &arg1 , &arg2 , &arg3 ),
                         new core::operation( c_kernel, &arg0 , &arg1 , &arg2 , &arg3 ));
}
