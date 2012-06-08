template< typename T0 , void(*fp)( T0& )>
struct kernel_inliner<void(*)( T0& ), fp> {
    static inline void latch(cfunctor* o, T0& arg0 ){
        info<T0>::typed::modify<0>(arg0, o);
    }
    static inline void invoke(sfunctor* o){
        fp( info<T0>::typed::revised<0>(o) );
    }
    static inline void cleanup(sfunctor* o){
        info<T0>::typed::deallocate<0>(o);
    }
    static inline void weight(cfunctor* o){
        info<T0>::typed::weight<0>(o);
    }
    static inline void place(sfunctor* o){
        info<T0>::typed::place<0>(o);
    }
};
template< typename T0 , typename T1 , void(*fp)( T0& , T1& )>
struct kernel_inliner<void(*)( T0& , T1& ), fp> {
    static inline void latch(cfunctor* o, T0& arg0 , T1& arg1 ){
        info<T0>::typed::modify<0>(arg0, o); info<T1>::typed::modify<1>(arg1, o);
    }
    static inline void invoke(sfunctor* o){
        fp( info<T0>::typed::revised<0>(o) , info<T1>::typed::revised<1>(o) );
    }
    static inline void cleanup(sfunctor* o){
        info<T0>::typed::deallocate<0>(o); info<T1>::typed::deallocate<1>(o);
    }
    static inline void weight(cfunctor* o){
        info<T0>::typed::weight<0>(o); info<T1>::typed::weight<1>(o);
    }
    static inline void place(sfunctor* o){
        info<T0>::typed::place<0>(o); info<T1>::typed::place<1>(o);
    }
};




template< typename T0 , typename T1 , typename T2 , void(*fp)( T0& , T1& , T2& )>
struct kernel_inliner<void(*)( T0& , T1& , T2& ), fp> {

    static inline void latch(cfunctor* o, T0& arg0 , T1& arg1 , T2& arg2 ){
        info<T0>::typed::modify<0>(arg0, o); info<T1>::typed::modify<1>(arg1, o); info<T2>::typed::modify<2>(arg2, o);
    }
    static inline void invoke(sfunctor* o){
        fp( info<T0>::typed::revised<0>(o) , info<T1>::typed::revised<1>(o) , info<T2>::typed::revised<2>(o) );
    }
    static inline void cleanup(sfunctor* o){
        info<T0>::typed::deallocate<0>(o); info<T1>::typed::deallocate<1>(o); info<T2>::typed::deallocate<2>(o);
    }
    static inline void weight(cfunctor* o){
        info<T0>::typed::weight<0>(o); info<T1>::typed::weight<1>(o); info<T2>::typed::weight<2>(o);
    }
    static inline void place(sfunctor* o){
        info<T0>::typed::place<0>(o); info<T1>::typed::place<1>(o); info<T2>::typed::place<2>(o);
    }
};
template< typename T0 , typename T1 , typename T2 , typename T3 , void(*fp)( T0& , T1& , T2& , T3& )>
struct kernel_inliner<void(*)( T0& , T1& , T2& , T3& ), fp> {
    static inline void latch(cfunctor* o, T0& arg0 , T1& arg1 , T2& arg2 , T3& arg3 ){
        info<T0>::typed::modify<0>(arg0, o); info<T1>::typed::modify<1>(arg1, o); info<T2>::typed::modify<2>(arg2, o); info<T3>::typed::modify<3>(arg3, o);
    }
    static inline void invoke(sfunctor* o){
        fp( info<T0>::typed::revised<0>(o) , info<T1>::typed::revised<1>(o) , info<T2>::typed::revised<2>(o) , info<T3>::typed::revised<3>(o) );
    }
    static inline void cleanup(sfunctor* o){
        info<T0>::typed::deallocate<0>(o); info<T1>::typed::deallocate<1>(o); info<T2>::typed::deallocate<2>(o); info<T3>::typed::deallocate<3>(o);
    }
    static inline void weight(cfunctor* o){
        info<T0>::typed::weight<0>(o); info<T1>::typed::weight<1>(o); info<T2>::typed::weight<2>(o); info<T3>::typed::weight<3>(o);
    }
    static inline void place(sfunctor* o){
        info<T0>::typed::place<0>(o); info<T1>::typed::place<1>(o); info<T2>::typed::place<2>(o); info<T3>::typed::place<3>(o);
    }
};
template< typename T0 , typename T1 , typename T2 , typename T3 , typename T4 , void(*fp)( T0& , T1& , T2& , T3& , T4& )>
struct kernel_inliner<void(*)( T0& , T1& , T2& , T3& , T4& ), fp> {
    static inline void latch(cfunctor* o, T0& arg0 , T1& arg1 , T2& arg2 , T3& arg3 , T4& arg4 ){
        info<T0>::typed::modify<0>(arg0, o); info<T1>::typed::modify<1>(arg1, o); info<T2>::typed::modify<2>(arg2, o); info<T3>::typed::modify<3>(arg3, o); info<T4>::typed::modify<4>(arg4, o);
    }
    static inline void invoke(sfunctor* o){
        fp( info<T0>::typed::revised<0>(o) , info<T1>::typed::revised<1>(o) , info<T2>::typed::revised<2>(o) , info<T3>::typed::revised<3>(o) , info<T4>::typed::revised<4>(o) );
    }
    static inline void cleanup(sfunctor* o){
        info<T0>::typed::deallocate<0>(o); info<T1>::typed::deallocate<1>(o); info<T2>::typed::deallocate<2>(o); info<T3>::typed::deallocate<3>(o); info<T4>::typed::deallocate<4>(o);
    }
    static inline void weight(cfunctor* o){
        info<T0>::typed::weight<0>(o); info<T1>::typed::weight<1>(o); info<T2>::typed::weight<2>(o); info<T3>::typed::weight<3>(o); info<T4>::typed::weight<4>(o);
    }
    static inline void place(sfunctor* o){
        info<T0>::typed::place<0>(o); info<T1>::typed::place<1>(o); info<T2>::typed::place<2>(o); info<T3>::typed::place<3>(o); info<T4>::typed::place<4>(o);
    }
};
template< typename T0 , typename T1 , typename T2 , typename T3 , typename T4 , typename T5 , void(*fp)( T0& , T1& , T2& , T3& , T4& , T5& )>
struct kernel_inliner<void(*)( T0& , T1& , T2& , T3& , T4& , T5& ), fp> {
    static inline void latch(cfunctor* o, T0& arg0 , T1& arg1 , T2& arg2 , T3& arg3 , T4& arg4 , T5& arg5 ){
        info<T0>::typed::modify<0>(arg0, o); info<T1>::typed::modify<1>(arg1, o); info<T2>::typed::modify<2>(arg2, o); info<T3>::typed::modify<3>(arg3, o); info<T4>::typed::modify<4>(arg4, o); info<T5>::typed::modify<5>(arg5, o);
    }
    static inline void invoke(sfunctor* o){
        fp( info<T0>::typed::revised<0>(o) , info<T1>::typed::revised<1>(o) , info<T2>::typed::revised<2>(o) , info<T3>::typed::revised<3>(o) , info<T4>::typed::revised<4>(o) , info<T5>::typed::revised<5>(o) );
    }
    static inline void cleanup(sfunctor* o){
        info<T0>::typed::deallocate<0>(o); info<T1>::typed::deallocate<1>(o); info<T2>::typed::deallocate<2>(o); info<T3>::typed::deallocate<3>(o); info<T4>::typed::deallocate<4>(o); info<T5>::typed::deallocate<5>(o);
    }
    static inline void weight(cfunctor* o){
        info<T0>::typed::weight<0>(o); info<T1>::typed::weight<1>(o); info<T2>::typed::weight<2>(o); info<T3>::typed::weight<3>(o); info<T4>::typed::weight<4>(o); info<T5>::typed::weight<5>(o);
    }
    static inline void place(sfunctor* o){
        info<T0>::typed::place<0>(o); info<T1>::typed::place<1>(o); info<T2>::typed::place<2>(o); info<T3>::typed::place<3>(o); info<T4>::typed::place<4>(o); info<T5>::typed::place<5>(o);
    }
};
template< typename T0 , typename T1 , typename T2 , typename T3 , typename T4 , typename T5 , typename T6 , void(*fp)( T0& , T1& , T2& , T3& , T4& , T5& , T6& )>
struct kernel_inliner<void(*)( T0& , T1& , T2& , T3& , T4& , T5& , T6& ), fp> {
    static inline void latch(cfunctor* o, T0& arg0 , T1& arg1 , T2& arg2 , T3& arg3 , T4& arg4 , T5& arg5 , T6& arg6 ){
        info<T0>::typed::modify<0>(arg0, o); info<T1>::typed::modify<1>(arg1, o); info<T2>::typed::modify<2>(arg2, o); info<T3>::typed::modify<3>(arg3, o); info<T4>::typed::modify<4>(arg4, o); info<T5>::typed::modify<5>(arg5, o); info<T6>::typed::modify<6>(arg6, o);
    }
    static inline void invoke(sfunctor* o){
        fp( info<T0>::typed::revised<0>(o) , info<T1>::typed::revised<1>(o) , info<T2>::typed::revised<2>(o) , info<T3>::typed::revised<3>(o) , info<T4>::typed::revised<4>(o) , info<T5>::typed::revised<5>(o) , info<T6>::typed::revised<6>(o) );
    }
    static inline void cleanup(sfunctor* o){
        info<T0>::typed::deallocate<0>(o); info<T1>::typed::deallocate<1>(o); info<T2>::typed::deallocate<2>(o); info<T3>::typed::deallocate<3>(o); info<T4>::typed::deallocate<4>(o); info<T5>::typed::deallocate<5>(o); info<T6>::typed::deallocate<6>(o);
    }
    static inline void weight(cfunctor* o){
        info<T0>::typed::weight<0>(o); info<T1>::typed::weight<1>(o); info<T2>::typed::weight<2>(o); info<T3>::typed::weight<3>(o); info<T4>::typed::weight<4>(o); info<T5>::typed::weight<5>(o); info<T6>::typed::weight<6>(o);
    }
    static inline void place(sfunctor* o){
        info<T0>::typed::place<0>(o); info<T1>::typed::place<1>(o); info<T2>::typed::place<2>(o); info<T3>::typed::place<3>(o); info<T4>::typed::place<4>(o); info<T5>::typed::place<5>(o); info<T6>::typed::place<6>(o);
    }
};
template< typename T0 , typename T1 , typename T2 , typename T3 , typename T4 , typename T5 , typename T6 , typename T7 , void(*fp)( T0& , T1& , T2& , T3& , T4& , T5& , T6& , T7& )>
struct kernel_inliner<void(*)( T0& , T1& , T2& , T3& , T4& , T5& , T6& , T7& ), fp> {
    static inline void latch(cfunctor* o, T0& arg0 , T1& arg1 , T2& arg2 , T3& arg3 , T4& arg4 , T5& arg5 , T6& arg6 , T7& arg7 ){
        info<T0>::typed::modify<0>(arg0, o); info<T1>::typed::modify<1>(arg1, o); info<T2>::typed::modify<2>(arg2, o); info<T3>::typed::modify<3>(arg3, o); info<T4>::typed::modify<4>(arg4, o); info<T5>::typed::modify<5>(arg5, o); info<T6>::typed::modify<6>(arg6, o); info<T7>::typed::modify<7>(arg7, o);
    }
    static inline void invoke(sfunctor* o){
        fp( info<T0>::typed::revised<0>(o) , info<T1>::typed::revised<1>(o) , info<T2>::typed::revised<2>(o) , info<T3>::typed::revised<3>(o) , info<T4>::typed::revised<4>(o) , info<T5>::typed::revised<5>(o) , info<T6>::typed::revised<6>(o) , info<T7>::typed::revised<7>(o) );
    }
    static inline void cleanup(sfunctor* o){
        info<T0>::typed::deallocate<0>(o); info<T1>::typed::deallocate<1>(o); info<T2>::typed::deallocate<2>(o); info<T3>::typed::deallocate<3>(o); info<T4>::typed::deallocate<4>(o); info<T5>::typed::deallocate<5>(o); info<T6>::typed::deallocate<6>(o); info<T7>::typed::deallocate<7>(o);
    }
    static inline void weight(cfunctor* o){
        info<T0>::typed::weight<0>(o); info<T1>::typed::weight<1>(o); info<T2>::typed::weight<2>(o); info<T3>::typed::weight<3>(o); info<T4>::typed::weight<4>(o); info<T5>::typed::weight<5>(o); info<T6>::typed::weight<6>(o); info<T7>::typed::weight<7>(o);
    }
    static inline void place(sfunctor* o){
        info<T0>::typed::place<0>(o); info<T1>::typed::place<1>(o); info<T2>::typed::place<2>(o); info<T3>::typed::place<3>(o); info<T4>::typed::place<4>(o); info<T5>::typed::place<5>(o); info<T6>::typed::place<6>(o); info<T7>::typed::place<7>(o);
    }
};
template< typename T0 , typename T1 , typename T2 , typename T3 , typename T4 , typename T5 , typename T6 , typename T7 , typename T8 , void(*fp)( T0& , T1& , T2& , T3& , T4& , T5& , T6& , T7& , T8& )>
struct kernel_inliner<void(*)( T0& , T1& , T2& , T3& , T4& , T5& , T6& , T7& , T8& ), fp> {
    static inline void latch(cfunctor* o, T0& arg0 , T1& arg1 , T2& arg2 , T3& arg3 , T4& arg4 , T5& arg5 , T6& arg6 , T7& arg7 , T8& arg8 ){
        info<T0>::typed::modify<0>(arg0, o); info<T1>::typed::modify<1>(arg1, o); info<T2>::typed::modify<2>(arg2, o); info<T3>::typed::modify<3>(arg3, o); info<T4>::typed::modify<4>(arg4, o); info<T5>::typed::modify<5>(arg5, o); info<T6>::typed::modify<6>(arg6, o); info<T7>::typed::modify<7>(arg7, o); info<T8>::typed::modify<8>(arg8, o);
    }
    static inline void invoke(sfunctor* o){
        fp( info<T0>::typed::revised<0>(o) , info<T1>::typed::revised<1>(o) , info<T2>::typed::revised<2>(o) , info<T3>::typed::revised<3>(o) , info<T4>::typed::revised<4>(o) , info<T5>::typed::revised<5>(o) , info<T6>::typed::revised<6>(o) , info<T7>::typed::revised<7>(o) , info<T8>::typed::revised<8>(o) );
    }
    static inline void cleanup(sfunctor* o){
        info<T0>::typed::deallocate<0>(o); info<T1>::typed::deallocate<1>(o); info<T2>::typed::deallocate<2>(o); info<T3>::typed::deallocate<3>(o); info<T4>::typed::deallocate<4>(o); info<T5>::typed::deallocate<5>(o); info<T6>::typed::deallocate<6>(o); info<T7>::typed::deallocate<7>(o); info<T8>::typed::deallocate<8>(o);
    }
    static inline void weight(cfunctor* o){
        info<T0>::typed::weight<0>(o); info<T1>::typed::weight<1>(o); info<T2>::typed::weight<2>(o); info<T3>::typed::weight<3>(o); info<T4>::typed::weight<4>(o); info<T5>::typed::weight<5>(o); info<T6>::typed::weight<6>(o); info<T7>::typed::weight<7>(o); info<T8>::typed::weight<8>(o);
    }
    static inline void place(sfunctor* o){
        info<T0>::typed::place<0>(o); info<T1>::typed::place<1>(o); info<T2>::typed::place<2>(o); info<T3>::typed::place<3>(o); info<T4>::typed::place<4>(o); info<T5>::typed::place<5>(o); info<T6>::typed::place<6>(o); info<T7>::typed::place<7>(o); info<T8>::typed::place<8>(o);
    }
};
template< typename T0 , typename T1 , typename T2 , typename T3 , typename T4 , typename T5 , typename T6 , typename T7 , typename T8 , typename T9 , void(*fp)( T0& , T1& , T2& , T3& , T4& , T5& , T6& , T7& , T8& , T9& )>
struct kernel_inliner<void(*)( T0& , T1& , T2& , T3& , T4& , T5& , T6& , T7& , T8& , T9& ), fp> {
    static inline void latch(cfunctor* o, T0& arg0 , T1& arg1 , T2& arg2 , T3& arg3 , T4& arg4 , T5& arg5 , T6& arg6 , T7& arg7 , T8& arg8 , T9& arg9 ){
        info<T0>::typed::modify<0>(arg0, o); info<T1>::typed::modify<1>(arg1, o); info<T2>::typed::modify<2>(arg2, o); info<T3>::typed::modify<3>(arg3, o); info<T4>::typed::modify<4>(arg4, o); info<T5>::typed::modify<5>(arg5, o); info<T6>::typed::modify<6>(arg6, o); info<T7>::typed::modify<7>(arg7, o); info<T8>::typed::modify<8>(arg8, o); info<T9>::typed::modify<9>(arg9, o);
    }
    static inline void invoke(sfunctor* o){
        fp( info<T0>::typed::revised<0>(o) , info<T1>::typed::revised<1>(o) , info<T2>::typed::revised<2>(o) , info<T3>::typed::revised<3>(o) , info<T4>::typed::revised<4>(o) , info<T5>::typed::revised<5>(o) , info<T6>::typed::revised<6>(o) , info<T7>::typed::revised<7>(o) , info<T8>::typed::revised<8>(o) , info<T9>::typed::revised<9>(o) );
    }
    static inline void cleanup(sfunctor* o){
        info<T0>::typed::deallocate<0>(o); info<T1>::typed::deallocate<1>(o); info<T2>::typed::deallocate<2>(o); info<T3>::typed::deallocate<3>(o); info<T4>::typed::deallocate<4>(o); info<T5>::typed::deallocate<5>(o); info<T6>::typed::deallocate<6>(o); info<T7>::typed::deallocate<7>(o); info<T8>::typed::deallocate<8>(o); info<T9>::typed::deallocate<9>(o);
    }
    static inline void weight(cfunctor* o){
        info<T0>::typed::weight<0>(o); info<T1>::typed::weight<1>(o); info<T2>::typed::weight<2>(o); info<T3>::typed::weight<3>(o); info<T4>::typed::weight<4>(o); info<T5>::typed::weight<5>(o); info<T6>::typed::weight<6>(o); info<T7>::typed::weight<7>(o); info<T8>::typed::weight<8>(o); info<T9>::typed::weight<9>(o);
    }
    static inline void place(sfunctor* o){
        info<T0>::typed::place<0>(o); info<T1>::typed::place<1>(o); info<T2>::typed::place<2>(o); info<T3>::typed::place<3>(o); info<T4>::typed::place<4>(o); info<T5>::typed::place<5>(o); info<T6>::typed::place<6>(o); info<T7>::typed::place<7>(o); info<T8>::typed::place<8>(o); info<T9>::typed::place<9>(o);
    }
};
