template<typename T0 , void(*fp)( T0& )>
struct kernel_inliner<void(*)( T0& ), fp> {
    static inline void latch(cfunctor* o, T0& arg0 ){
        info<T0>::typed::modify<0>(arg0, o);
        if(info<T0>::typed::pin<0>(o)) return;
        ambient::controller.submit(o);
    }
    static inline void invoke(sfunctor* o){
        (*fp)( info<T0>::typed::revised<0>(o) );
    }
    static inline void cleanup(sfunctor* o){
        info<T0>::typed::deallocate<0>(o);
    }
    static inline void deploy(sfunctor* o, size_t target){
        info<T0>::typed::deploy<0>(o,target);
    }
    static inline bool ready(sfunctor* o){
        return (info<T0>::typed::ready<0>(o) && true);
    }
};
template<typename T0 , typename T1 , void(*fp)( T0& , T1& )>
struct kernel_inliner<void(*)( T0& , T1& ), fp> {
    static inline void latch(cfunctor* o, T0& arg0 , T1& arg1 ){
        info<T0>::typed::modify<0>(arg0, o); info<T1>::typed::modify<1>(arg1, o);
        if(info<T0>::typed::pin<0>(o)) return; if(info<T1>::typed::pin<1>(o)) return;
        ambient::controller.submit(o);
    }
    static inline void invoke(sfunctor* o){
        (*fp)( info<T0>::typed::revised<0>(o) , info<T1>::typed::revised<1>(o) );
    }
    static inline void cleanup(sfunctor* o){
        info<T0>::typed::deallocate<0>(o); info<T1>::typed::deallocate<1>(o);
    }
    static inline void deploy(sfunctor* o, size_t target){
        info<T0>::typed::deploy<0>(o,target); info<T1>::typed::deploy<1>(o,target);
    }
    static inline bool ready(sfunctor* o){
        return (info<T0>::typed::ready<0>(o) && info<T1>::typed::ready<1>(o) && true);
    }
};



template<typename T0 , typename T1 , typename T2 , void(*fp)( T0& , T1& , T2& )>
struct kernel_inliner<void(*)( T0& , T1& , T2& ), fp> {

    static inline void latch(cfunctor* o, T0& arg0 , T1& arg1 , T2& arg2 ){
        info<T0>::typed::modify<0>(arg0, o); info<T1>::typed::modify<1>(arg1, o); info<T2>::typed::modify<2>(arg2, o);
        if(info<T0>::typed::pin<0>(o)) return; if(info<T1>::typed::pin<1>(o)) return; if(info<T2>::typed::pin<2>(o)) return;
        ambient::controller.submit(o);
    }
    static inline void invoke(sfunctor* o){
        (*fp)( info<T0>::typed::revised<0>(o) , info<T1>::typed::revised<1>(o) , info<T2>::typed::revised<2>(o) );
    }
    static inline void cleanup(sfunctor* o){
        info<T0>::typed::deallocate<0>(o); info<T1>::typed::deallocate<1>(o); info<T2>::typed::deallocate<2>(o);
    }
    static inline void deploy(sfunctor* o, size_t target){
        info<T0>::typed::deploy<0>(o,target); info<T1>::typed::deploy<1>(o,target); info<T2>::typed::deploy<2>(o,target);
    }
    static inline bool ready(sfunctor* o){
        return (info<T0>::typed::ready<0>(o) && info<T1>::typed::ready<1>(o) && info<T2>::typed::ready<2>(o) && true);
    }
};
template<typename T0 , typename T1 , typename T2 , typename T3 , void(*fp)( T0& , T1& , T2& , T3& )>
struct kernel_inliner<void(*)( T0& , T1& , T2& , T3& ), fp> {
    static inline void latch(cfunctor* o, T0& arg0 , T1& arg1 , T2& arg2 , T3& arg3 ){
        info<T0>::typed::modify<0>(arg0, o); info<T1>::typed::modify<1>(arg1, o); info<T2>::typed::modify<2>(arg2, o); info<T3>::typed::modify<3>(arg3, o);
        if(info<T0>::typed::pin<0>(o)) return; if(info<T1>::typed::pin<1>(o)) return; if(info<T2>::typed::pin<2>(o)) return; if(info<T3>::typed::pin<3>(o)) return;
        ambient::controller.submit(o);
    }
    static inline void invoke(sfunctor* o){
        (*fp)( info<T0>::typed::revised<0>(o) , info<T1>::typed::revised<1>(o) , info<T2>::typed::revised<2>(o) , info<T3>::typed::revised<3>(o) );
    }
    static inline void cleanup(sfunctor* o){
        info<T0>::typed::deallocate<0>(o); info<T1>::typed::deallocate<1>(o); info<T2>::typed::deallocate<2>(o); info<T3>::typed::deallocate<3>(o);
    }
    static inline void deploy(sfunctor* o, size_t target){
        info<T0>::typed::deploy<0>(o,target); info<T1>::typed::deploy<1>(o,target); info<T2>::typed::deploy<2>(o,target); info<T3>::typed::deploy<3>(o,target);
    }
    static inline bool ready(sfunctor* o){
        return (info<T0>::typed::ready<0>(o) && info<T1>::typed::ready<1>(o) && info<T2>::typed::ready<2>(o) && info<T3>::typed::ready<3>(o) && true);
    }
};
template<typename T0 , typename T1 , typename T2 , typename T3 , typename T4 , void(*fp)( T0& , T1& , T2& , T3& , T4& )>
struct kernel_inliner<void(*)( T0& , T1& , T2& , T3& , T4& ), fp> {
    static inline void latch(cfunctor* o, T0& arg0 , T1& arg1 , T2& arg2 , T3& arg3 , T4& arg4 ){
        info<T0>::typed::modify<0>(arg0, o); info<T1>::typed::modify<1>(arg1, o); info<T2>::typed::modify<2>(arg2, o); info<T3>::typed::modify<3>(arg3, o); info<T4>::typed::modify<4>(arg4, o);
        if(info<T0>::typed::pin<0>(o)) return; if(info<T1>::typed::pin<1>(o)) return; if(info<T2>::typed::pin<2>(o)) return; if(info<T3>::typed::pin<3>(o)) return; if(info<T4>::typed::pin<4>(o)) return;
        ambient::controller.submit(o);
    }
    static inline void invoke(sfunctor* o){
        (*fp)( info<T0>::typed::revised<0>(o) , info<T1>::typed::revised<1>(o) , info<T2>::typed::revised<2>(o) , info<T3>::typed::revised<3>(o) , info<T4>::typed::revised<4>(o) );
    }
    static inline void cleanup(sfunctor* o){
        info<T0>::typed::deallocate<0>(o); info<T1>::typed::deallocate<1>(o); info<T2>::typed::deallocate<2>(o); info<T3>::typed::deallocate<3>(o); info<T4>::typed::deallocate<4>(o);
    }
    static inline void deploy(sfunctor* o, size_t target){
        info<T0>::typed::deploy<0>(o,target); info<T1>::typed::deploy<1>(o,target); info<T2>::typed::deploy<2>(o,target); info<T3>::typed::deploy<3>(o,target); info<T4>::typed::deploy<4>(o,target);
    }
    static inline bool ready(sfunctor* o){
        return (info<T0>::typed::ready<0>(o) && info<T1>::typed::ready<1>(o) && info<T2>::typed::ready<2>(o) && info<T3>::typed::ready<3>(o) && info<T4>::typed::ready<4>(o) && true);
    }
};
template<typename T0 , typename T1 , typename T2 , typename T3 , typename T4 , typename T5 , void(*fp)( T0& , T1& , T2& , T3& , T4& , T5& )>
struct kernel_inliner<void(*)( T0& , T1& , T2& , T3& , T4& , T5& ), fp> {
    static inline void latch(cfunctor* o, T0& arg0 , T1& arg1 , T2& arg2 , T3& arg3 , T4& arg4 , T5& arg5 ){
        info<T0>::typed::modify<0>(arg0, o); info<T1>::typed::modify<1>(arg1, o); info<T2>::typed::modify<2>(arg2, o); info<T3>::typed::modify<3>(arg3, o); info<T4>::typed::modify<4>(arg4, o); info<T5>::typed::modify<5>(arg5, o);
        if(info<T0>::typed::pin<0>(o)) return; if(info<T1>::typed::pin<1>(o)) return; if(info<T2>::typed::pin<2>(o)) return; if(info<T3>::typed::pin<3>(o)) return; if(info<T4>::typed::pin<4>(o)) return; if(info<T5>::typed::pin<5>(o)) return;
        ambient::controller.submit(o);
    }
    static inline void invoke(sfunctor* o){
        (*fp)( info<T0>::typed::revised<0>(o) , info<T1>::typed::revised<1>(o) , info<T2>::typed::revised<2>(o) , info<T3>::typed::revised<3>(o) , info<T4>::typed::revised<4>(o) , info<T5>::typed::revised<5>(o) );
    }
    static inline void cleanup(sfunctor* o){
        info<T0>::typed::deallocate<0>(o); info<T1>::typed::deallocate<1>(o); info<T2>::typed::deallocate<2>(o); info<T3>::typed::deallocate<3>(o); info<T4>::typed::deallocate<4>(o); info<T5>::typed::deallocate<5>(o);
    }
    static inline void deploy(sfunctor* o, size_t target){
        info<T0>::typed::deploy<0>(o,target); info<T1>::typed::deploy<1>(o,target); info<T2>::typed::deploy<2>(o,target); info<T3>::typed::deploy<3>(o,target); info<T4>::typed::deploy<4>(o,target); info<T5>::typed::deploy<5>(o,target);
    }
    static inline bool ready(sfunctor* o){
        return (info<T0>::typed::ready<0>(o) && info<T1>::typed::ready<1>(o) && info<T2>::typed::ready<2>(o) && info<T3>::typed::ready<3>(o) && info<T4>::typed::ready<4>(o) && info<T5>::typed::ready<5>(o) && true);
    }
};
template<typename T0 , typename T1 , typename T2 , typename T3 , typename T4 , typename T5 , typename T6 , void(*fp)( T0& , T1& , T2& , T3& , T4& , T5& , T6& )>
struct kernel_inliner<void(*)( T0& , T1& , T2& , T3& , T4& , T5& , T6& ), fp> {
    static inline void latch(cfunctor* o, T0& arg0 , T1& arg1 , T2& arg2 , T3& arg3 , T4& arg4 , T5& arg5 , T6& arg6 ){
        info<T0>::typed::modify<0>(arg0, o); info<T1>::typed::modify<1>(arg1, o); info<T2>::typed::modify<2>(arg2, o); info<T3>::typed::modify<3>(arg3, o); info<T4>::typed::modify<4>(arg4, o); info<T5>::typed::modify<5>(arg5, o); info<T6>::typed::modify<6>(arg6, o);
        if(info<T0>::typed::pin<0>(o)) return; if(info<T1>::typed::pin<1>(o)) return; if(info<T2>::typed::pin<2>(o)) return; if(info<T3>::typed::pin<3>(o)) return; if(info<T4>::typed::pin<4>(o)) return; if(info<T5>::typed::pin<5>(o)) return; if(info<T6>::typed::pin<6>(o)) return;
        ambient::controller.submit(o);
    }
    static inline void invoke(sfunctor* o){
        (*fp)( info<T0>::typed::revised<0>(o) , info<T1>::typed::revised<1>(o) , info<T2>::typed::revised<2>(o) , info<T3>::typed::revised<3>(o) , info<T4>::typed::revised<4>(o) , info<T5>::typed::revised<5>(o) , info<T6>::typed::revised<6>(o) );
    }
    static inline void cleanup(sfunctor* o){
        info<T0>::typed::deallocate<0>(o); info<T1>::typed::deallocate<1>(o); info<T2>::typed::deallocate<2>(o); info<T3>::typed::deallocate<3>(o); info<T4>::typed::deallocate<4>(o); info<T5>::typed::deallocate<5>(o); info<T6>::typed::deallocate<6>(o);
    }
    static inline void deploy(sfunctor* o, size_t target){
        info<T0>::typed::deploy<0>(o,target); info<T1>::typed::deploy<1>(o,target); info<T2>::typed::deploy<2>(o,target); info<T3>::typed::deploy<3>(o,target); info<T4>::typed::deploy<4>(o,target); info<T5>::typed::deploy<5>(o,target); info<T6>::typed::deploy<6>(o,target);
    }
    static inline bool ready(sfunctor* o){
        return (info<T0>::typed::ready<0>(o) && info<T1>::typed::ready<1>(o) && info<T2>::typed::ready<2>(o) && info<T3>::typed::ready<3>(o) && info<T4>::typed::ready<4>(o) && info<T5>::typed::ready<5>(o) && info<T6>::typed::ready<6>(o) && true);
    }
};
template<typename T0 , typename T1 , typename T2 , typename T3 , typename T4 , typename T5 , typename T6 , typename T7 , void(*fp)( T0& , T1& , T2& , T3& , T4& , T5& , T6& , T7& )>
struct kernel_inliner<void(*)( T0& , T1& , T2& , T3& , T4& , T5& , T6& , T7& ), fp> {
    static inline void latch(cfunctor* o, T0& arg0 , T1& arg1 , T2& arg2 , T3& arg3 , T4& arg4 , T5& arg5 , T6& arg6 , T7& arg7 ){
        info<T0>::typed::modify<0>(arg0, o); info<T1>::typed::modify<1>(arg1, o); info<T2>::typed::modify<2>(arg2, o); info<T3>::typed::modify<3>(arg3, o); info<T4>::typed::modify<4>(arg4, o); info<T5>::typed::modify<5>(arg5, o); info<T6>::typed::modify<6>(arg6, o); info<T7>::typed::modify<7>(arg7, o);
        if(info<T0>::typed::pin<0>(o)) return; if(info<T1>::typed::pin<1>(o)) return; if(info<T2>::typed::pin<2>(o)) return; if(info<T3>::typed::pin<3>(o)) return; if(info<T4>::typed::pin<4>(o)) return; if(info<T5>::typed::pin<5>(o)) return; if(info<T6>::typed::pin<6>(o)) return; if(info<T7>::typed::pin<7>(o)) return;
        ambient::controller.submit(o);
    }
    static inline void invoke(sfunctor* o){
        (*fp)( info<T0>::typed::revised<0>(o) , info<T1>::typed::revised<1>(o) , info<T2>::typed::revised<2>(o) , info<T3>::typed::revised<3>(o) , info<T4>::typed::revised<4>(o) , info<T5>::typed::revised<5>(o) , info<T6>::typed::revised<6>(o) , info<T7>::typed::revised<7>(o) );
    }
    static inline void cleanup(sfunctor* o){
        info<T0>::typed::deallocate<0>(o); info<T1>::typed::deallocate<1>(o); info<T2>::typed::deallocate<2>(o); info<T3>::typed::deallocate<3>(o); info<T4>::typed::deallocate<4>(o); info<T5>::typed::deallocate<5>(o); info<T6>::typed::deallocate<6>(o); info<T7>::typed::deallocate<7>(o);
    }
    static inline void deploy(sfunctor* o, size_t target){
        info<T0>::typed::deploy<0>(o,target); info<T1>::typed::deploy<1>(o,target); info<T2>::typed::deploy<2>(o,target); info<T3>::typed::deploy<3>(o,target); info<T4>::typed::deploy<4>(o,target); info<T5>::typed::deploy<5>(o,target); info<T6>::typed::deploy<6>(o,target); info<T7>::typed::deploy<7>(o,target);
    }
    static inline bool ready(sfunctor* o){
        return (info<T0>::typed::ready<0>(o) && info<T1>::typed::ready<1>(o) && info<T2>::typed::ready<2>(o) && info<T3>::typed::ready<3>(o) && info<T4>::typed::ready<4>(o) && info<T5>::typed::ready<5>(o) && info<T6>::typed::ready<6>(o) && info<T7>::typed::ready<7>(o) && true);
    }
};
template<typename T0 , typename T1 , typename T2 , typename T3 , typename T4 , typename T5 , typename T6 , typename T7 , typename T8 , void(*fp)( T0& , T1& , T2& , T3& , T4& , T5& , T6& , T7& , T8& )>
struct kernel_inliner<void(*)( T0& , T1& , T2& , T3& , T4& , T5& , T6& , T7& , T8& ), fp> {
    static inline void latch(cfunctor* o, T0& arg0 , T1& arg1 , T2& arg2 , T3& arg3 , T4& arg4 , T5& arg5 , T6& arg6 , T7& arg7 , T8& arg8 ){
        info<T0>::typed::modify<0>(arg0, o); info<T1>::typed::modify<1>(arg1, o); info<T2>::typed::modify<2>(arg2, o); info<T3>::typed::modify<3>(arg3, o); info<T4>::typed::modify<4>(arg4, o); info<T5>::typed::modify<5>(arg5, o); info<T6>::typed::modify<6>(arg6, o); info<T7>::typed::modify<7>(arg7, o); info<T8>::typed::modify<8>(arg8, o);
        if(info<T0>::typed::pin<0>(o)) return; if(info<T1>::typed::pin<1>(o)) return; if(info<T2>::typed::pin<2>(o)) return; if(info<T3>::typed::pin<3>(o)) return; if(info<T4>::typed::pin<4>(o)) return; if(info<T5>::typed::pin<5>(o)) return; if(info<T6>::typed::pin<6>(o)) return; if(info<T7>::typed::pin<7>(o)) return; if(info<T8>::typed::pin<8>(o)) return;
        ambient::controller.submit(o);
    }
    static inline void invoke(sfunctor* o){
        (*fp)( info<T0>::typed::revised<0>(o) , info<T1>::typed::revised<1>(o) , info<T2>::typed::revised<2>(o) , info<T3>::typed::revised<3>(o) , info<T4>::typed::revised<4>(o) , info<T5>::typed::revised<5>(o) , info<T6>::typed::revised<6>(o) , info<T7>::typed::revised<7>(o) , info<T8>::typed::revised<8>(o) );
    }
    static inline void cleanup(sfunctor* o){
        info<T0>::typed::deallocate<0>(o); info<T1>::typed::deallocate<1>(o); info<T2>::typed::deallocate<2>(o); info<T3>::typed::deallocate<3>(o); info<T4>::typed::deallocate<4>(o); info<T5>::typed::deallocate<5>(o); info<T6>::typed::deallocate<6>(o); info<T7>::typed::deallocate<7>(o); info<T8>::typed::deallocate<8>(o);
    }
    static inline void deploy(sfunctor* o, size_t target){
        info<T0>::typed::deploy<0>(o,target); info<T1>::typed::deploy<1>(o,target); info<T2>::typed::deploy<2>(o,target); info<T3>::typed::deploy<3>(o,target); info<T4>::typed::deploy<4>(o,target); info<T5>::typed::deploy<5>(o,target); info<T6>::typed::deploy<6>(o,target); info<T7>::typed::deploy<7>(o,target); info<T8>::typed::deploy<8>(o,target);
    }
    static inline bool ready(sfunctor* o){
        return (info<T0>::typed::ready<0>(o) && info<T1>::typed::ready<1>(o) && info<T2>::typed::ready<2>(o) && info<T3>::typed::ready<3>(o) && info<T4>::typed::ready<4>(o) && info<T5>::typed::ready<5>(o) && info<T6>::typed::ready<6>(o) && info<T7>::typed::ready<7>(o) && info<T8>::typed::ready<8>(o) && true);
    }
};
template<typename T0 , typename T1 , typename T2 , typename T3 , typename T4 , typename T5 , typename T6 , typename T7 , typename T8 , typename T9 , void(*fp)( T0& , T1& , T2& , T3& , T4& , T5& , T6& , T7& , T8& , T9& )>
struct kernel_inliner<void(*)( T0& , T1& , T2& , T3& , T4& , T5& , T6& , T7& , T8& , T9& ), fp> {
    static inline void latch(cfunctor* o, T0& arg0 , T1& arg1 , T2& arg2 , T3& arg3 , T4& arg4 , T5& arg5 , T6& arg6 , T7& arg7 , T8& arg8 , T9& arg9 ){
        info<T0>::typed::modify<0>(arg0, o); info<T1>::typed::modify<1>(arg1, o); info<T2>::typed::modify<2>(arg2, o); info<T3>::typed::modify<3>(arg3, o); info<T4>::typed::modify<4>(arg4, o); info<T5>::typed::modify<5>(arg5, o); info<T6>::typed::modify<6>(arg6, o); info<T7>::typed::modify<7>(arg7, o); info<T8>::typed::modify<8>(arg8, o); info<T9>::typed::modify<9>(arg9, o);
        if(info<T0>::typed::pin<0>(o)) return; if(info<T1>::typed::pin<1>(o)) return; if(info<T2>::typed::pin<2>(o)) return; if(info<T3>::typed::pin<3>(o)) return; if(info<T4>::typed::pin<4>(o)) return; if(info<T5>::typed::pin<5>(o)) return; if(info<T6>::typed::pin<6>(o)) return; if(info<T7>::typed::pin<7>(o)) return; if(info<T8>::typed::pin<8>(o)) return; if(info<T9>::typed::pin<9>(o)) return;
        ambient::controller.submit(o);
    }
    static inline void invoke(sfunctor* o){
        (*fp)( info<T0>::typed::revised<0>(o) , info<T1>::typed::revised<1>(o) , info<T2>::typed::revised<2>(o) , info<T3>::typed::revised<3>(o) , info<T4>::typed::revised<4>(o) , info<T5>::typed::revised<5>(o) , info<T6>::typed::revised<6>(o) , info<T7>::typed::revised<7>(o) , info<T8>::typed::revised<8>(o) , info<T9>::typed::revised<9>(o) );
    }
    static inline void cleanup(sfunctor* o){
        info<T0>::typed::deallocate<0>(o); info<T1>::typed::deallocate<1>(o); info<T2>::typed::deallocate<2>(o); info<T3>::typed::deallocate<3>(o); info<T4>::typed::deallocate<4>(o); info<T5>::typed::deallocate<5>(o); info<T6>::typed::deallocate<6>(o); info<T7>::typed::deallocate<7>(o); info<T8>::typed::deallocate<8>(o); info<T9>::typed::deallocate<9>(o);
    }
    static inline void deploy(sfunctor* o, size_t target){
        info<T0>::typed::deploy<0>(o,target); info<T1>::typed::deploy<1>(o,target); info<T2>::typed::deploy<2>(o,target); info<T3>::typed::deploy<3>(o,target); info<T4>::typed::deploy<4>(o,target); info<T5>::typed::deploy<5>(o,target); info<T6>::typed::deploy<6>(o,target); info<T7>::typed::deploy<7>(o,target); info<T8>::typed::deploy<8>(o,target); info<T9>::typed::deploy<9>(o,target);
    }
    static inline bool ready(sfunctor* o){
        return (info<T0>::typed::ready<0>(o) && info<T1>::typed::ready<1>(o) && info<T2>::typed::ready<2>(o) && info<T3>::typed::ready<3>(o) && info<T4>::typed::ready<4>(o) && info<T5>::typed::ready<5>(o) && info<T6>::typed::ready<6>(o) && info<T7>::typed::ready<7>(o) && info<T8>::typed::ready<8>(o) && info<T9>::typed::ready<9>(o) && true);
    }
};
template<typename T0 , typename T1 , typename T2 , typename T3 , typename T4 , typename T5 , typename T6 , typename T7 , typename T8 , typename T9 , typename T10 , void(*fp)( T0& , T1& , T2& , T3& , T4& , T5& , T6& , T7& , T8& , T9& , T10& )>
struct kernel_inliner<void(*)( T0& , T1& , T2& , T3& , T4& , T5& , T6& , T7& , T8& , T9& , T10& ), fp> {
    static inline void latch(cfunctor* o, T0& arg0 , T1& arg1 , T2& arg2 , T3& arg3 , T4& arg4 , T5& arg5 , T6& arg6 , T7& arg7 , T8& arg8 , T9& arg9 , T10& arg10 ){
        info<T0>::typed::modify<0>(arg0, o); info<T1>::typed::modify<1>(arg1, o); info<T2>::typed::modify<2>(arg2, o); info<T3>::typed::modify<3>(arg3, o); info<T4>::typed::modify<4>(arg4, o); info<T5>::typed::modify<5>(arg5, o); info<T6>::typed::modify<6>(arg6, o); info<T7>::typed::modify<7>(arg7, o); info<T8>::typed::modify<8>(arg8, o); info<T9>::typed::modify<9>(arg9, o); info<T10>::typed::modify<10>(arg10, o);
        if(info<T0>::typed::pin<0>(o)) return; if(info<T1>::typed::pin<1>(o)) return; if(info<T2>::typed::pin<2>(o)) return; if(info<T3>::typed::pin<3>(o)) return; if(info<T4>::typed::pin<4>(o)) return; if(info<T5>::typed::pin<5>(o)) return; if(info<T6>::typed::pin<6>(o)) return; if(info<T7>::typed::pin<7>(o)) return; if(info<T8>::typed::pin<8>(o)) return; if(info<T9>::typed::pin<9>(o)) return; if(info<T10>::typed::pin<10>(o)) return;
        ambient::controller.submit(o);
    }
    static inline void invoke(sfunctor* o){
        (*fp)( info<T0>::typed::revised<0>(o) , info<T1>::typed::revised<1>(o) , info<T2>::typed::revised<2>(o) , info<T3>::typed::revised<3>(o) , info<T4>::typed::revised<4>(o) , info<T5>::typed::revised<5>(o) , info<T6>::typed::revised<6>(o) , info<T7>::typed::revised<7>(o) , info<T8>::typed::revised<8>(o) , info<T9>::typed::revised<9>(o) , info<T10>::typed::revised<10>(o) );
    }
    static inline void cleanup(sfunctor* o){
        info<T0>::typed::deallocate<0>(o); info<T1>::typed::deallocate<1>(o); info<T2>::typed::deallocate<2>(o); info<T3>::typed::deallocate<3>(o); info<T4>::typed::deallocate<4>(o); info<T5>::typed::deallocate<5>(o); info<T6>::typed::deallocate<6>(o); info<T7>::typed::deallocate<7>(o); info<T8>::typed::deallocate<8>(o); info<T9>::typed::deallocate<9>(o); info<T10>::typed::deallocate<10>(o);
    }
    static inline void deploy(sfunctor* o, size_t target){
        info<T0>::typed::deploy<0>(o,target); info<T1>::typed::deploy<1>(o,target); info<T2>::typed::deploy<2>(o,target); info<T3>::typed::deploy<3>(o,target); info<T4>::typed::deploy<4>(o,target); info<T5>::typed::deploy<5>(o,target); info<T6>::typed::deploy<6>(o,target); info<T7>::typed::deploy<7>(o,target); info<T8>::typed::deploy<8>(o,target); info<T9>::typed::deploy<9>(o,target); info<T10>::typed::deploy<10>(o,target);
    }
    static inline bool ready(sfunctor* o){
        return (info<T0>::typed::ready<0>(o) && info<T1>::typed::ready<1>(o) && info<T2>::typed::ready<2>(o) && info<T3>::typed::ready<3>(o) && info<T4>::typed::ready<4>(o) && info<T5>::typed::ready<5>(o) && info<T6>::typed::ready<6>(o) && info<T7>::typed::ready<7>(o) && info<T8>::typed::ready<8>(o) && info<T9>::typed::ready<9>(o) && info<T10>::typed::ready<10>(o) && true);
    }
};
