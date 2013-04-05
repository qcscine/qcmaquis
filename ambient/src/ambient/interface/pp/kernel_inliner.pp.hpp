template< typename T0 , void(*fp)( T0& )>
struct kernel_inliner<void(*)( T0& ), fp> {
    typedef T0 t0;
    static const int arity = 1;
    template<complexity O>
    static inline void latch(cfunctor* o, T0& arg0 ){
        if(ambient::controller.tunable()){
            info<T0>::typed::template score<0>(arg0);
            ambient::controller.schedule<O>();
        }
        if(ambient::controller.remote()){
            info<T0>::typed::template modify_remote<0>(arg0);
            return;
        }else if(ambient::controller.local()){
            info<T0>::typed::template modify_local<0>(arg0, o);
        }else{
            info<T0>::typed::template modify<0>(arg0, o);
        }
        info<T0>::typed::template pin<0>(o) ||
        ambient::controller.queue(o);
    }
    static inline void invoke(cfunctor* o){
        (*fp)( info<T0>::typed::template revised<0>(o) );
    }
    static inline void cleanup(cfunctor* o){
        info<T0>::typed::template deallocate<0>(o);
    }
    static inline bool ready(cfunctor* o){
        return (info<T0>::typed::template ready<0>(o) && true);
    }
};
template< typename T0 , typename T1 , void(*fp)( T0& , T1& )>
struct kernel_inliner<void(*)( T0& , T1& ), fp> {
    typedef T0 t0; typedef T1 t1;
    static const int arity = 2;
    template<complexity O>
    static inline void latch(cfunctor* o, T0& arg0 , T1& arg1 ){
        if(ambient::controller.tunable()){
            info<T0>::typed::template score<0>(arg0); info<T1>::typed::template score<1>(arg1);
            ambient::controller.schedule<O>();
        }
        if(ambient::controller.remote()){
            info<T0>::typed::template modify_remote<0>(arg0); info<T1>::typed::template modify_remote<1>(arg1);
            return;
        }else if(ambient::controller.local()){
            info<T0>::typed::template modify_local<0>(arg0, o); info<T1>::typed::template modify_local<1>(arg1, o);
        }else{
            info<T0>::typed::template modify<0>(arg0, o); info<T1>::typed::template modify<1>(arg1, o);
        }
        info<T0>::typed::template pin<0>(o) || info<T1>::typed::template pin<1>(o) ||
        ambient::controller.queue(o);
    }
    static inline void invoke(cfunctor* o){
        (*fp)( info<T0>::typed::template revised<0>(o) , info<T1>::typed::template revised<1>(o) );
    }
    static inline void cleanup(cfunctor* o){
        info<T0>::typed::template deallocate<0>(o); info<T1>::typed::template deallocate<1>(o);
    }
    static inline bool ready(cfunctor* o){
        return (info<T0>::typed::template ready<0>(o) && info<T1>::typed::template ready<1>(o) && true);
    }
};
template< typename T0 , typename T1 , typename T2 , void(*fp)( T0& , T1& , T2& )>
struct kernel_inliner<void(*)( T0& , T1& , T2& ), fp> {
    typedef T0 t0; typedef T1 t1; typedef T2 t2;
    static const int arity = 3;
    template<complexity O>
    static inline void latch(cfunctor* o, T0& arg0 , T1& arg1 , T2& arg2 ){
        if(ambient::controller.tunable()){
            info<T0>::typed::template score<0>(arg0); info<T1>::typed::template score<1>(arg1); info<T2>::typed::template score<2>(arg2);
            ambient::controller.schedule<O>();
        }
        if(ambient::controller.remote()){
            info<T0>::typed::template modify_remote<0>(arg0); info<T1>::typed::template modify_remote<1>(arg1); info<T2>::typed::template modify_remote<2>(arg2);
            return;
        }else if(ambient::controller.local()){
            info<T0>::typed::template modify_local<0>(arg0, o); info<T1>::typed::template modify_local<1>(arg1, o); info<T2>::typed::template modify_local<2>(arg2, o);
        }else{
            info<T0>::typed::template modify<0>(arg0, o); info<T1>::typed::template modify<1>(arg1, o); info<T2>::typed::template modify<2>(arg2, o);
        }
        info<T0>::typed::template pin<0>(o) || info<T1>::typed::template pin<1>(o) || info<T2>::typed::template pin<2>(o) ||
        ambient::controller.queue(o);
    }
    static inline void invoke(cfunctor* o){
        (*fp)( info<T0>::typed::template revised<0>(o) , info<T1>::typed::template revised<1>(o) , info<T2>::typed::template revised<2>(o) );
    }
    static inline void cleanup(cfunctor* o){
        info<T0>::typed::template deallocate<0>(o); info<T1>::typed::template deallocate<1>(o); info<T2>::typed::template deallocate<2>(o);
    }
    static inline bool ready(cfunctor* o){
        return (info<T0>::typed::template ready<0>(o) && info<T1>::typed::template ready<1>(o) && info<T2>::typed::template ready<2>(o) && true);
    }
};
template< typename T0 , typename T1 , typename T2 , typename T3 , void(*fp)( T0& , T1& , T2& , T3& )>
struct kernel_inliner<void(*)( T0& , T1& , T2& , T3& ), fp> {
    typedef T0 t0; typedef T1 t1; typedef T2 t2; typedef T3 t3;
    static const int arity = 4;
    template<complexity O>
    static inline void latch(cfunctor* o, T0& arg0 , T1& arg1 , T2& arg2 , T3& arg3 ){
        if(ambient::controller.tunable()){
            info<T0>::typed::template score<0>(arg0); info<T1>::typed::template score<1>(arg1); info<T2>::typed::template score<2>(arg2); info<T3>::typed::template score<3>(arg3);
            ambient::controller.schedule<O>();
        }
        if(ambient::controller.remote()){
            info<T0>::typed::template modify_remote<0>(arg0); info<T1>::typed::template modify_remote<1>(arg1); info<T2>::typed::template modify_remote<2>(arg2); info<T3>::typed::template modify_remote<3>(arg3);
            return;
        }else if(ambient::controller.local()){
            info<T0>::typed::template modify_local<0>(arg0, o); info<T1>::typed::template modify_local<1>(arg1, o); info<T2>::typed::template modify_local<2>(arg2, o); info<T3>::typed::template modify_local<3>(arg3, o);
        }else{
            info<T0>::typed::template modify<0>(arg0, o); info<T1>::typed::template modify<1>(arg1, o); info<T2>::typed::template modify<2>(arg2, o); info<T3>::typed::template modify<3>(arg3, o);
        }
        info<T0>::typed::template pin<0>(o) || info<T1>::typed::template pin<1>(o) || info<T2>::typed::template pin<2>(o) || info<T3>::typed::template pin<3>(o) ||
        ambient::controller.queue(o);
    }
    static inline void invoke(cfunctor* o){
        (*fp)( info<T0>::typed::template revised<0>(o) , info<T1>::typed::template revised<1>(o) , info<T2>::typed::template revised<2>(o) , info<T3>::typed::template revised<3>(o) );
    }
    static inline void cleanup(cfunctor* o){
        info<T0>::typed::template deallocate<0>(o); info<T1>::typed::template deallocate<1>(o); info<T2>::typed::template deallocate<2>(o); info<T3>::typed::template deallocate<3>(o);
    }
    static inline bool ready(cfunctor* o){
        return (info<T0>::typed::template ready<0>(o) && info<T1>::typed::template ready<1>(o) && info<T2>::typed::template ready<2>(o) && info<T3>::typed::template ready<3>(o) && true);
    }
};
template< typename T0 , typename T1 , typename T2 , typename T3 , typename T4 , void(*fp)( T0& , T1& , T2& , T3& , T4& )>
struct kernel_inliner<void(*)( T0& , T1& , T2& , T3& , T4& ), fp> {
    typedef T0 t0; typedef T1 t1; typedef T2 t2; typedef T3 t3; typedef T4 t4;
    static const int arity = 5;
    template<complexity O>
    static inline void latch(cfunctor* o, T0& arg0 , T1& arg1 , T2& arg2 , T3& arg3 , T4& arg4 ){
        if(ambient::controller.tunable()){
            info<T0>::typed::template score<0>(arg0); info<T1>::typed::template score<1>(arg1); info<T2>::typed::template score<2>(arg2); info<T3>::typed::template score<3>(arg3); info<T4>::typed::template score<4>(arg4);
            ambient::controller.schedule<O>();
        }
        if(ambient::controller.remote()){
            info<T0>::typed::template modify_remote<0>(arg0); info<T1>::typed::template modify_remote<1>(arg1); info<T2>::typed::template modify_remote<2>(arg2); info<T3>::typed::template modify_remote<3>(arg3); info<T4>::typed::template modify_remote<4>(arg4);
            return;
        }else if(ambient::controller.local()){
            info<T0>::typed::template modify_local<0>(arg0, o); info<T1>::typed::template modify_local<1>(arg1, o); info<T2>::typed::template modify_local<2>(arg2, o); info<T3>::typed::template modify_local<3>(arg3, o); info<T4>::typed::template modify_local<4>(arg4, o);
        }else{
            info<T0>::typed::template modify<0>(arg0, o); info<T1>::typed::template modify<1>(arg1, o); info<T2>::typed::template modify<2>(arg2, o); info<T3>::typed::template modify<3>(arg3, o); info<T4>::typed::template modify<4>(arg4, o);
        }
        info<T0>::typed::template pin<0>(o) || info<T1>::typed::template pin<1>(o) || info<T2>::typed::template pin<2>(o) || info<T3>::typed::template pin<3>(o) || info<T4>::typed::template pin<4>(o) ||
        ambient::controller.queue(o);
    }
    static inline void invoke(cfunctor* o){
        (*fp)( info<T0>::typed::template revised<0>(o) , info<T1>::typed::template revised<1>(o) , info<T2>::typed::template revised<2>(o) , info<T3>::typed::template revised<3>(o) , info<T4>::typed::template revised<4>(o) );
    }
    static inline void cleanup(cfunctor* o){
        info<T0>::typed::template deallocate<0>(o); info<T1>::typed::template deallocate<1>(o); info<T2>::typed::template deallocate<2>(o); info<T3>::typed::template deallocate<3>(o); info<T4>::typed::template deallocate<4>(o);
    }
    static inline bool ready(cfunctor* o){
        return (info<T0>::typed::template ready<0>(o) && info<T1>::typed::template ready<1>(o) && info<T2>::typed::template ready<2>(o) && info<T3>::typed::template ready<3>(o) && info<T4>::typed::template ready<4>(o) && true);
    }
};
template< typename T0 , typename T1 , typename T2 , typename T3 , typename T4 , typename T5 , void(*fp)( T0& , T1& , T2& , T3& , T4& , T5& )>
struct kernel_inliner<void(*)( T0& , T1& , T2& , T3& , T4& , T5& ), fp> {
    typedef T0 t0; typedef T1 t1; typedef T2 t2; typedef T3 t3; typedef T4 t4; typedef T5 t5;
    static const int arity = 6;
    template<complexity O>
    static inline void latch(cfunctor* o, T0& arg0 , T1& arg1 , T2& arg2 , T3& arg3 , T4& arg4 , T5& arg5 ){
        if(ambient::controller.tunable()){
            info<T0>::typed::template score<0>(arg0); info<T1>::typed::template score<1>(arg1); info<T2>::typed::template score<2>(arg2); info<T3>::typed::template score<3>(arg3); info<T4>::typed::template score<4>(arg4); info<T5>::typed::template score<5>(arg5);
            ambient::controller.schedule<O>();
        }
        if(ambient::controller.remote()){
            info<T0>::typed::template modify_remote<0>(arg0); info<T1>::typed::template modify_remote<1>(arg1); info<T2>::typed::template modify_remote<2>(arg2); info<T3>::typed::template modify_remote<3>(arg3); info<T4>::typed::template modify_remote<4>(arg4); info<T5>::typed::template modify_remote<5>(arg5);
            return;
        }else if(ambient::controller.local()){
            info<T0>::typed::template modify_local<0>(arg0, o); info<T1>::typed::template modify_local<1>(arg1, o); info<T2>::typed::template modify_local<2>(arg2, o); info<T3>::typed::template modify_local<3>(arg3, o); info<T4>::typed::template modify_local<4>(arg4, o); info<T5>::typed::template modify_local<5>(arg5, o);
        }else{
            info<T0>::typed::template modify<0>(arg0, o); info<T1>::typed::template modify<1>(arg1, o); info<T2>::typed::template modify<2>(arg2, o); info<T3>::typed::template modify<3>(arg3, o); info<T4>::typed::template modify<4>(arg4, o); info<T5>::typed::template modify<5>(arg5, o);
        }
        info<T0>::typed::template pin<0>(o) || info<T1>::typed::template pin<1>(o) || info<T2>::typed::template pin<2>(o) || info<T3>::typed::template pin<3>(o) || info<T4>::typed::template pin<4>(o) || info<T5>::typed::template pin<5>(o) ||
        ambient::controller.queue(o);
    }
    static inline void invoke(cfunctor* o){
        (*fp)( info<T0>::typed::template revised<0>(o) , info<T1>::typed::template revised<1>(o) , info<T2>::typed::template revised<2>(o) , info<T3>::typed::template revised<3>(o) , info<T4>::typed::template revised<4>(o) , info<T5>::typed::template revised<5>(o) );
    }
    static inline void cleanup(cfunctor* o){
        info<T0>::typed::template deallocate<0>(o); info<T1>::typed::template deallocate<1>(o); info<T2>::typed::template deallocate<2>(o); info<T3>::typed::template deallocate<3>(o); info<T4>::typed::template deallocate<4>(o); info<T5>::typed::template deallocate<5>(o);
    }
    static inline bool ready(cfunctor* o){
        return (info<T0>::typed::template ready<0>(o) && info<T1>::typed::template ready<1>(o) && info<T2>::typed::template ready<2>(o) && info<T3>::typed::template ready<3>(o) && info<T4>::typed::template ready<4>(o) && info<T5>::typed::template ready<5>(o) && true);
    }
};
template< typename T0 , typename T1 , typename T2 , typename T3 , typename T4 , typename T5 , typename T6 , void(*fp)( T0& , T1& , T2& , T3& , T4& , T5& , T6& )>
struct kernel_inliner<void(*)( T0& , T1& , T2& , T3& , T4& , T5& , T6& ), fp> {
    typedef T0 t0; typedef T1 t1; typedef T2 t2; typedef T3 t3; typedef T4 t4; typedef T5 t5; typedef T6 t6;
    static const int arity = 7;
    template<complexity O>
    static inline void latch(cfunctor* o, T0& arg0 , T1& arg1 , T2& arg2 , T3& arg3 , T4& arg4 , T5& arg5 , T6& arg6 ){
        if(ambient::controller.tunable()){
            info<T0>::typed::template score<0>(arg0); info<T1>::typed::template score<1>(arg1); info<T2>::typed::template score<2>(arg2); info<T3>::typed::template score<3>(arg3); info<T4>::typed::template score<4>(arg4); info<T5>::typed::template score<5>(arg5); info<T6>::typed::template score<6>(arg6);
            ambient::controller.schedule<O>();
        }
        if(ambient::controller.remote()){
            info<T0>::typed::template modify_remote<0>(arg0); info<T1>::typed::template modify_remote<1>(arg1); info<T2>::typed::template modify_remote<2>(arg2); info<T3>::typed::template modify_remote<3>(arg3); info<T4>::typed::template modify_remote<4>(arg4); info<T5>::typed::template modify_remote<5>(arg5); info<T6>::typed::template modify_remote<6>(arg6);
            return;
        }else if(ambient::controller.local()){
            info<T0>::typed::template modify_local<0>(arg0, o); info<T1>::typed::template modify_local<1>(arg1, o); info<T2>::typed::template modify_local<2>(arg2, o); info<T3>::typed::template modify_local<3>(arg3, o); info<T4>::typed::template modify_local<4>(arg4, o); info<T5>::typed::template modify_local<5>(arg5, o); info<T6>::typed::template modify_local<6>(arg6, o);
        }else{
            info<T0>::typed::template modify<0>(arg0, o); info<T1>::typed::template modify<1>(arg1, o); info<T2>::typed::template modify<2>(arg2, o); info<T3>::typed::template modify<3>(arg3, o); info<T4>::typed::template modify<4>(arg4, o); info<T5>::typed::template modify<5>(arg5, o); info<T6>::typed::template modify<6>(arg6, o);
        }
        info<T0>::typed::template pin<0>(o) || info<T1>::typed::template pin<1>(o) || info<T2>::typed::template pin<2>(o) || info<T3>::typed::template pin<3>(o) || info<T4>::typed::template pin<4>(o) || info<T5>::typed::template pin<5>(o) || info<T6>::typed::template pin<6>(o) ||
        ambient::controller.queue(o);
    }
    static inline void invoke(cfunctor* o){
        (*fp)( info<T0>::typed::template revised<0>(o) , info<T1>::typed::template revised<1>(o) , info<T2>::typed::template revised<2>(o) , info<T3>::typed::template revised<3>(o) , info<T4>::typed::template revised<4>(o) , info<T5>::typed::template revised<5>(o) , info<T6>::typed::template revised<6>(o) );
    }
    static inline void cleanup(cfunctor* o){
        info<T0>::typed::template deallocate<0>(o); info<T1>::typed::template deallocate<1>(o); info<T2>::typed::template deallocate<2>(o); info<T3>::typed::template deallocate<3>(o); info<T4>::typed::template deallocate<4>(o); info<T5>::typed::template deallocate<5>(o); info<T6>::typed::template deallocate<6>(o);
    }
    static inline bool ready(cfunctor* o){
        return (info<T0>::typed::template ready<0>(o) && info<T1>::typed::template ready<1>(o) && info<T2>::typed::template ready<2>(o) && info<T3>::typed::template ready<3>(o) && info<T4>::typed::template ready<4>(o) && info<T5>::typed::template ready<5>(o) && info<T6>::typed::template ready<6>(o) && true);
    }
};
template< typename T0 , typename T1 , typename T2 , typename T3 , typename T4 , typename T5 , typename T6 , typename T7 , void(*fp)( T0& , T1& , T2& , T3& , T4& , T5& , T6& , T7& )>
struct kernel_inliner<void(*)( T0& , T1& , T2& , T3& , T4& , T5& , T6& , T7& ), fp> {
    typedef T0 t0; typedef T1 t1; typedef T2 t2; typedef T3 t3; typedef T4 t4; typedef T5 t5; typedef T6 t6; typedef T7 t7;
    static const int arity = 8;
    template<complexity O>
    static inline void latch(cfunctor* o, T0& arg0 , T1& arg1 , T2& arg2 , T3& arg3 , T4& arg4 , T5& arg5 , T6& arg6 , T7& arg7 ){
        if(ambient::controller.tunable()){
            info<T0>::typed::template score<0>(arg0); info<T1>::typed::template score<1>(arg1); info<T2>::typed::template score<2>(arg2); info<T3>::typed::template score<3>(arg3); info<T4>::typed::template score<4>(arg4); info<T5>::typed::template score<5>(arg5); info<T6>::typed::template score<6>(arg6); info<T7>::typed::template score<7>(arg7);
            ambient::controller.schedule<O>();
        }
        if(ambient::controller.remote()){
            info<T0>::typed::template modify_remote<0>(arg0); info<T1>::typed::template modify_remote<1>(arg1); info<T2>::typed::template modify_remote<2>(arg2); info<T3>::typed::template modify_remote<3>(arg3); info<T4>::typed::template modify_remote<4>(arg4); info<T5>::typed::template modify_remote<5>(arg5); info<T6>::typed::template modify_remote<6>(arg6); info<T7>::typed::template modify_remote<7>(arg7);
            return;
        }else if(ambient::controller.local()){
            info<T0>::typed::template modify_local<0>(arg0, o); info<T1>::typed::template modify_local<1>(arg1, o); info<T2>::typed::template modify_local<2>(arg2, o); info<T3>::typed::template modify_local<3>(arg3, o); info<T4>::typed::template modify_local<4>(arg4, o); info<T5>::typed::template modify_local<5>(arg5, o); info<T6>::typed::template modify_local<6>(arg6, o); info<T7>::typed::template modify_local<7>(arg7, o);
        }else{
            info<T0>::typed::template modify<0>(arg0, o); info<T1>::typed::template modify<1>(arg1, o); info<T2>::typed::template modify<2>(arg2, o); info<T3>::typed::template modify<3>(arg3, o); info<T4>::typed::template modify<4>(arg4, o); info<T5>::typed::template modify<5>(arg5, o); info<T6>::typed::template modify<6>(arg6, o); info<T7>::typed::template modify<7>(arg7, o);
        }
        info<T0>::typed::template pin<0>(o) || info<T1>::typed::template pin<1>(o) || info<T2>::typed::template pin<2>(o) || info<T3>::typed::template pin<3>(o) || info<T4>::typed::template pin<4>(o) || info<T5>::typed::template pin<5>(o) || info<T6>::typed::template pin<6>(o) || info<T7>::typed::template pin<7>(o) ||
        ambient::controller.queue(o);
    }
    static inline void invoke(cfunctor* o){
        (*fp)( info<T0>::typed::template revised<0>(o) , info<T1>::typed::template revised<1>(o) , info<T2>::typed::template revised<2>(o) , info<T3>::typed::template revised<3>(o) , info<T4>::typed::template revised<4>(o) , info<T5>::typed::template revised<5>(o) , info<T6>::typed::template revised<6>(o) , info<T7>::typed::template revised<7>(o) );
    }
    static inline void cleanup(cfunctor* o){
        info<T0>::typed::template deallocate<0>(o); info<T1>::typed::template deallocate<1>(o); info<T2>::typed::template deallocate<2>(o); info<T3>::typed::template deallocate<3>(o); info<T4>::typed::template deallocate<4>(o); info<T5>::typed::template deallocate<5>(o); info<T6>::typed::template deallocate<6>(o); info<T7>::typed::template deallocate<7>(o);
    }
    static inline bool ready(cfunctor* o){
        return (info<T0>::typed::template ready<0>(o) && info<T1>::typed::template ready<1>(o) && info<T2>::typed::template ready<2>(o) && info<T3>::typed::template ready<3>(o) && info<T4>::typed::template ready<4>(o) && info<T5>::typed::template ready<5>(o) && info<T6>::typed::template ready<6>(o) && info<T7>::typed::template ready<7>(o) && true);
    }
};
template< typename T0 , typename T1 , typename T2 , typename T3 , typename T4 , typename T5 , typename T6 , typename T7 , typename T8 , void(*fp)( T0& , T1& , T2& , T3& , T4& , T5& , T6& , T7& , T8& )>
struct kernel_inliner<void(*)( T0& , T1& , T2& , T3& , T4& , T5& , T6& , T7& , T8& ), fp> {
    typedef T0 t0; typedef T1 t1; typedef T2 t2; typedef T3 t3; typedef T4 t4; typedef T5 t5; typedef T6 t6; typedef T7 t7; typedef T8 t8;
    static const int arity = 9;
    template<complexity O>
    static inline void latch(cfunctor* o, T0& arg0 , T1& arg1 , T2& arg2 , T3& arg3 , T4& arg4 , T5& arg5 , T6& arg6 , T7& arg7 , T8& arg8 ){
        if(ambient::controller.tunable()){
            info<T0>::typed::template score<0>(arg0); info<T1>::typed::template score<1>(arg1); info<T2>::typed::template score<2>(arg2); info<T3>::typed::template score<3>(arg3); info<T4>::typed::template score<4>(arg4); info<T5>::typed::template score<5>(arg5); info<T6>::typed::template score<6>(arg6); info<T7>::typed::template score<7>(arg7); info<T8>::typed::template score<8>(arg8);
            ambient::controller.schedule<O>();
        }
        if(ambient::controller.remote()){
            info<T0>::typed::template modify_remote<0>(arg0); info<T1>::typed::template modify_remote<1>(arg1); info<T2>::typed::template modify_remote<2>(arg2); info<T3>::typed::template modify_remote<3>(arg3); info<T4>::typed::template modify_remote<4>(arg4); info<T5>::typed::template modify_remote<5>(arg5); info<T6>::typed::template modify_remote<6>(arg6); info<T7>::typed::template modify_remote<7>(arg7); info<T8>::typed::template modify_remote<8>(arg8);
            return;
        }else if(ambient::controller.local()){
            info<T0>::typed::template modify_local<0>(arg0, o); info<T1>::typed::template modify_local<1>(arg1, o); info<T2>::typed::template modify_local<2>(arg2, o); info<T3>::typed::template modify_local<3>(arg3, o); info<T4>::typed::template modify_local<4>(arg4, o); info<T5>::typed::template modify_local<5>(arg5, o); info<T6>::typed::template modify_local<6>(arg6, o); info<T7>::typed::template modify_local<7>(arg7, o); info<T8>::typed::template modify_local<8>(arg8, o);
        }else{
            info<T0>::typed::template modify<0>(arg0, o); info<T1>::typed::template modify<1>(arg1, o); info<T2>::typed::template modify<2>(arg2, o); info<T3>::typed::template modify<3>(arg3, o); info<T4>::typed::template modify<4>(arg4, o); info<T5>::typed::template modify<5>(arg5, o); info<T6>::typed::template modify<6>(arg6, o); info<T7>::typed::template modify<7>(arg7, o); info<T8>::typed::template modify<8>(arg8, o);
        }
        info<T0>::typed::template pin<0>(o) || info<T1>::typed::template pin<1>(o) || info<T2>::typed::template pin<2>(o) || info<T3>::typed::template pin<3>(o) || info<T4>::typed::template pin<4>(o) || info<T5>::typed::template pin<5>(o) || info<T6>::typed::template pin<6>(o) || info<T7>::typed::template pin<7>(o) || info<T8>::typed::template pin<8>(o) ||
        ambient::controller.queue(o);
    }
    static inline void invoke(cfunctor* o){
        (*fp)( info<T0>::typed::template revised<0>(o) , info<T1>::typed::template revised<1>(o) , info<T2>::typed::template revised<2>(o) , info<T3>::typed::template revised<3>(o) , info<T4>::typed::template revised<4>(o) , info<T5>::typed::template revised<5>(o) , info<T6>::typed::template revised<6>(o) , info<T7>::typed::template revised<7>(o) , info<T8>::typed::template revised<8>(o) );
    }
    static inline void cleanup(cfunctor* o){
        info<T0>::typed::template deallocate<0>(o); info<T1>::typed::template deallocate<1>(o); info<T2>::typed::template deallocate<2>(o); info<T3>::typed::template deallocate<3>(o); info<T4>::typed::template deallocate<4>(o); info<T5>::typed::template deallocate<5>(o); info<T6>::typed::template deallocate<6>(o); info<T7>::typed::template deallocate<7>(o); info<T8>::typed::template deallocate<8>(o);
    }
    static inline bool ready(cfunctor* o){
        return (info<T0>::typed::template ready<0>(o) && info<T1>::typed::template ready<1>(o) && info<T2>::typed::template ready<2>(o) && info<T3>::typed::template ready<3>(o) && info<T4>::typed::template ready<4>(o) && info<T5>::typed::template ready<5>(o) && info<T6>::typed::template ready<6>(o) && info<T7>::typed::template ready<7>(o) && info<T8>::typed::template ready<8>(o) && true);
    }
};
template< typename T0 , typename T1 , typename T2 , typename T3 , typename T4 , typename T5 , typename T6 , typename T7 , typename T8 , typename T9 , void(*fp)( T0& , T1& , T2& , T3& , T4& , T5& , T6& , T7& , T8& , T9& )>
struct kernel_inliner<void(*)( T0& , T1& , T2& , T3& , T4& , T5& , T6& , T7& , T8& , T9& ), fp> {
    typedef T0 t0; typedef T1 t1; typedef T2 t2; typedef T3 t3; typedef T4 t4; typedef T5 t5; typedef T6 t6; typedef T7 t7; typedef T8 t8; typedef T9 t9;
    static const int arity = 10;
    template<complexity O>
    static inline void latch(cfunctor* o, T0& arg0 , T1& arg1 , T2& arg2 , T3& arg3 , T4& arg4 , T5& arg5 , T6& arg6 , T7& arg7 , T8& arg8 , T9& arg9 ){
        if(ambient::controller.tunable()){
            info<T0>::typed::template score<0>(arg0); info<T1>::typed::template score<1>(arg1); info<T2>::typed::template score<2>(arg2); info<T3>::typed::template score<3>(arg3); info<T4>::typed::template score<4>(arg4); info<T5>::typed::template score<5>(arg5); info<T6>::typed::template score<6>(arg6); info<T7>::typed::template score<7>(arg7); info<T8>::typed::template score<8>(arg8); info<T9>::typed::template score<9>(arg9);
            ambient::controller.schedule<O>();
        }
        if(ambient::controller.remote()){
            info<T0>::typed::template modify_remote<0>(arg0); info<T1>::typed::template modify_remote<1>(arg1); info<T2>::typed::template modify_remote<2>(arg2); info<T3>::typed::template modify_remote<3>(arg3); info<T4>::typed::template modify_remote<4>(arg4); info<T5>::typed::template modify_remote<5>(arg5); info<T6>::typed::template modify_remote<6>(arg6); info<T7>::typed::template modify_remote<7>(arg7); info<T8>::typed::template modify_remote<8>(arg8); info<T9>::typed::template modify_remote<9>(arg9);
            return;
        }else if(ambient::controller.local()){
            info<T0>::typed::template modify_local<0>(arg0, o); info<T1>::typed::template modify_local<1>(arg1, o); info<T2>::typed::template modify_local<2>(arg2, o); info<T3>::typed::template modify_local<3>(arg3, o); info<T4>::typed::template modify_local<4>(arg4, o); info<T5>::typed::template modify_local<5>(arg5, o); info<T6>::typed::template modify_local<6>(arg6, o); info<T7>::typed::template modify_local<7>(arg7, o); info<T8>::typed::template modify_local<8>(arg8, o); info<T9>::typed::template modify_local<9>(arg9, o);
        }else{
            info<T0>::typed::template modify<0>(arg0, o); info<T1>::typed::template modify<1>(arg1, o); info<T2>::typed::template modify<2>(arg2, o); info<T3>::typed::template modify<3>(arg3, o); info<T4>::typed::template modify<4>(arg4, o); info<T5>::typed::template modify<5>(arg5, o); info<T6>::typed::template modify<6>(arg6, o); info<T7>::typed::template modify<7>(arg7, o); info<T8>::typed::template modify<8>(arg8, o); info<T9>::typed::template modify<9>(arg9, o);
        }
        info<T0>::typed::template pin<0>(o) || info<T1>::typed::template pin<1>(o) || info<T2>::typed::template pin<2>(o) || info<T3>::typed::template pin<3>(o) || info<T4>::typed::template pin<4>(o) || info<T5>::typed::template pin<5>(o) || info<T6>::typed::template pin<6>(o) || info<T7>::typed::template pin<7>(o) || info<T8>::typed::template pin<8>(o) || info<T9>::typed::template pin<9>(o) ||
        ambient::controller.queue(o);
    }
    static inline void invoke(cfunctor* o){
        (*fp)( info<T0>::typed::template revised<0>(o) , info<T1>::typed::template revised<1>(o) , info<T2>::typed::template revised<2>(o) , info<T3>::typed::template revised<3>(o) , info<T4>::typed::template revised<4>(o) , info<T5>::typed::template revised<5>(o) , info<T6>::typed::template revised<6>(o) , info<T7>::typed::template revised<7>(o) , info<T8>::typed::template revised<8>(o) , info<T9>::typed::template revised<9>(o) );
    }
    static inline void cleanup(cfunctor* o){
        info<T0>::typed::template deallocate<0>(o); info<T1>::typed::template deallocate<1>(o); info<T2>::typed::template deallocate<2>(o); info<T3>::typed::template deallocate<3>(o); info<T4>::typed::template deallocate<4>(o); info<T5>::typed::template deallocate<5>(o); info<T6>::typed::template deallocate<6>(o); info<T7>::typed::template deallocate<7>(o); info<T8>::typed::template deallocate<8>(o); info<T9>::typed::template deallocate<9>(o);
    }
    static inline bool ready(cfunctor* o){
        return (info<T0>::typed::template ready<0>(o) && info<T1>::typed::template ready<1>(o) && info<T2>::typed::template ready<2>(o) && info<T3>::typed::template ready<3>(o) && info<T4>::typed::template ready<4>(o) && info<T5>::typed::template ready<5>(o) && info<T6>::typed::template ready<6>(o) && info<T7>::typed::template ready<7>(o) && info<T8>::typed::template ready<8>(o) && info<T9>::typed::template ready<9>(o) && true);
    }
};
template< typename T0 , typename T1 , typename T2 , typename T3 , typename T4 , typename T5 , typename T6 , typename T7 , typename T8 , typename T9 , typename T10 , void(*fp)( T0& , T1& , T2& , T3& , T4& , T5& , T6& , T7& , T8& , T9& , T10& )>
struct kernel_inliner<void(*)( T0& , T1& , T2& , T3& , T4& , T5& , T6& , T7& , T8& , T9& , T10& ), fp> {
    typedef T0 t0; typedef T1 t1; typedef T2 t2; typedef T3 t3; typedef T4 t4; typedef T5 t5; typedef T6 t6; typedef T7 t7; typedef T8 t8; typedef T9 t9; typedef T10 t10;
    static const int arity = 11;
    template<complexity O>
    static inline void latch(cfunctor* o, T0& arg0 , T1& arg1 , T2& arg2 , T3& arg3 , T4& arg4 , T5& arg5 , T6& arg6 , T7& arg7 , T8& arg8 , T9& arg9 , T10& arg10 ){
        if(ambient::controller.tunable()){
            info<T0>::typed::template score<0>(arg0); info<T1>::typed::template score<1>(arg1); info<T2>::typed::template score<2>(arg2); info<T3>::typed::template score<3>(arg3); info<T4>::typed::template score<4>(arg4); info<T5>::typed::template score<5>(arg5); info<T6>::typed::template score<6>(arg6); info<T7>::typed::template score<7>(arg7); info<T8>::typed::template score<8>(arg8); info<T9>::typed::template score<9>(arg9); info<T10>::typed::template score<10>(arg10);
            ambient::controller.schedule<O>();
        }
        if(ambient::controller.remote()){
            info<T0>::typed::template modify_remote<0>(arg0); info<T1>::typed::template modify_remote<1>(arg1); info<T2>::typed::template modify_remote<2>(arg2); info<T3>::typed::template modify_remote<3>(arg3); info<T4>::typed::template modify_remote<4>(arg4); info<T5>::typed::template modify_remote<5>(arg5); info<T6>::typed::template modify_remote<6>(arg6); info<T7>::typed::template modify_remote<7>(arg7); info<T8>::typed::template modify_remote<8>(arg8); info<T9>::typed::template modify_remote<9>(arg9); info<T10>::typed::template modify_remote<10>(arg10);
            return;
        }else if(ambient::controller.local()){
            info<T0>::typed::template modify_local<0>(arg0, o); info<T1>::typed::template modify_local<1>(arg1, o); info<T2>::typed::template modify_local<2>(arg2, o); info<T3>::typed::template modify_local<3>(arg3, o); info<T4>::typed::template modify_local<4>(arg4, o); info<T5>::typed::template modify_local<5>(arg5, o); info<T6>::typed::template modify_local<6>(arg6, o); info<T7>::typed::template modify_local<7>(arg7, o); info<T8>::typed::template modify_local<8>(arg8, o); info<T9>::typed::template modify_local<9>(arg9, o); info<T10>::typed::template modify_local<10>(arg10, o);
        }else{
            info<T0>::typed::template modify<0>(arg0, o); info<T1>::typed::template modify<1>(arg1, o); info<T2>::typed::template modify<2>(arg2, o); info<T3>::typed::template modify<3>(arg3, o); info<T4>::typed::template modify<4>(arg4, o); info<T5>::typed::template modify<5>(arg5, o); info<T6>::typed::template modify<6>(arg6, o); info<T7>::typed::template modify<7>(arg7, o); info<T8>::typed::template modify<8>(arg8, o); info<T9>::typed::template modify<9>(arg9, o); info<T10>::typed::template modify<10>(arg10, o);
        }
        info<T0>::typed::template pin<0>(o) || info<T1>::typed::template pin<1>(o) || info<T2>::typed::template pin<2>(o) || info<T3>::typed::template pin<3>(o) || info<T4>::typed::template pin<4>(o) || info<T5>::typed::template pin<5>(o) || info<T6>::typed::template pin<6>(o) || info<T7>::typed::template pin<7>(o) || info<T8>::typed::template pin<8>(o) || info<T9>::typed::template pin<9>(o) || info<T10>::typed::template pin<10>(o) ||
        ambient::controller.queue(o);
    }
    static inline void invoke(cfunctor* o){
        (*fp)( info<T0>::typed::template revised<0>(o) , info<T1>::typed::template revised<1>(o) , info<T2>::typed::template revised<2>(o) , info<T3>::typed::template revised<3>(o) , info<T4>::typed::template revised<4>(o) , info<T5>::typed::template revised<5>(o) , info<T6>::typed::template revised<6>(o) , info<T7>::typed::template revised<7>(o) , info<T8>::typed::template revised<8>(o) , info<T9>::typed::template revised<9>(o) , info<T10>::typed::template revised<10>(o) );
    }
    static inline void cleanup(cfunctor* o){
        info<T0>::typed::template deallocate<0>(o); info<T1>::typed::template deallocate<1>(o); info<T2>::typed::template deallocate<2>(o); info<T3>::typed::template deallocate<3>(o); info<T4>::typed::template deallocate<4>(o); info<T5>::typed::template deallocate<5>(o); info<T6>::typed::template deallocate<6>(o); info<T7>::typed::template deallocate<7>(o); info<T8>::typed::template deallocate<8>(o); info<T9>::typed::template deallocate<9>(o); info<T10>::typed::template deallocate<10>(o);
    }
    static inline bool ready(cfunctor* o){
        return (info<T0>::typed::template ready<0>(o) && info<T1>::typed::template ready<1>(o) && info<T2>::typed::template ready<2>(o) && info<T3>::typed::template ready<3>(o) && info<T4>::typed::template ready<4>(o) && info<T5>::typed::template ready<5>(o) && info<T6>::typed::template ready<6>(o) && info<T7>::typed::template ready<7>(o) && info<T8>::typed::template ready<8>(o) && info<T9>::typed::template ready<9>(o) && info<T10>::typed::template ready<10>(o) && true);
    }
};
