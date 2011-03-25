template< typename FP, typename T0 >
operation( FP op, T0 *arg0 );
template < typename T0 >
void extract_template(void (*)( T0& ));
template < typename T0 >
void prototype_template(void (*)( T0& ));
template < typename T0 > void extract_template(void (*)( pinned T0& )); template < typename T0 > void prototype_template(void (*)( pinned T0& ));
template< typename FP, typename T0 , typename T1 >
operation( FP op, T0 *arg0 , T1 *arg1 );
template < typename T0 , typename T1 >
void extract_template(void (*)( T0& , T1& ));
template < typename T0 , typename T1 >
void prototype_template(void (*)( T0& , T1& ));
template < typename T0 , typename T1 > void extract_template(void (*)( pinned T0& , T1& )); template < typename T0 , typename T1 > void prototype_template(void (*)( pinned T0& , T1& )); template < typename T0 , typename T1 > void extract_template(void (*)( T0& , pinned T1& )); template < typename T0 , typename T1 > void prototype_template(void (*)( T0& , pinned T1& ));
template< typename FP, typename T0 , typename T1 , typename T2 >
operation( FP op, T0 *arg0 , T1 *arg1 , T2 *arg2 );
template < typename T0 , typename T1 , typename T2 >
void extract_template(void (*)( T0& , T1& , T2& ));
template < typename T0 , typename T1 , typename T2 >
void prototype_template(void (*)( T0& , T1& , T2& ));
template < typename T0 , typename T1 , typename T2 > void extract_template(void (*)( pinned T0& , T1& , T2& )); template < typename T0 , typename T1 , typename T2 > void prototype_template(void (*)( pinned T0& , T1& , T2& )); template < typename T0 , typename T1 , typename T2 > void extract_template(void (*)( T0& , pinned T1& , T2& )); template < typename T0 , typename T1 , typename T2 > void prototype_template(void (*)( T0& , pinned T1& , T2& )); template < typename T0 , typename T1 , typename T2 > void extract_template(void (*)( T0& , T1& , pinned T2& )); template < typename T0 , typename T1 , typename T2 > void prototype_template(void (*)( T0& , T1& , pinned T2& ));







template< typename FP, typename T0 , typename T1 , typename T2 , typename T3 >
operation( FP op, T0 *arg0 , T1 *arg1 , T2 *arg2 , T3 *arg3 );
template < typename T0 , typename T1 , typename T2 , typename T3 >
void extract_template(void (*)( T0& , T1& , T2& , T3& ));
template < typename T0 , typename T1 , typename T2 , typename T3 >
void prototype_template(void (*)( T0& , T1& , T2& , T3& ));
template < typename T0 , typename T1 , typename T2 , typename T3 > void extract_template(void (*)( pinned T0& , T1& , T2& , T3& )); template < typename T0 , typename T1 , typename T2 , typename T3 > void prototype_template(void (*)( pinned T0& , T1& , T2& , T3& )); template < typename T0 , typename T1 , typename T2 , typename T3 > void extract_template(void (*)( T0& , pinned T1& , T2& , T3& )); template < typename T0 , typename T1 , typename T2 , typename T3 > void prototype_template(void (*)( T0& , pinned T1& , T2& , T3& )); template < typename T0 , typename T1 , typename T2 , typename T3 > void extract_template(void (*)( T0& , T1& , pinned T2& , T3& )); template < typename T0 , typename T1 , typename T2 , typename T3 > void prototype_template(void (*)( T0& , T1& , pinned T2& , T3& )); template < typename T0 , typename T1 , typename T2 , typename T3 > void extract_template(void (*)( T0& , T1& , T2& , pinned T3& )); template < typename T0 , typename T1 , typename T2 , typename T3 > void prototype_template(void (*)( T0& , T1& , T2& , pinned T3& ));
