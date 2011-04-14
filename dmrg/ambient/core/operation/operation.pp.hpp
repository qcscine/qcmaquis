namespace core {
template< typename FP, typename T0 >
operation::operation( FP op, T0 *arg0 ){
    this->init();
    this->operation_ptr = (void(*)())op;
    this->count = 1;
    this->profiles = (p_profile**)malloc(sizeof(p_profile*)*this->count);
    this->arguments = new boost::shared_ptr<void>[this->count];
    this->arguments[0] = boost::shared_ptr<void>(get_handle(*arg0));
    void(operation::*ptr)(FP); ptr = &operation::prototype_template;
    this->prototype = (void(operation::*)())ptr;
    ptr = &operation::extract_template;
    this->extract = (void(operation::*)())ptr;
}
template < typename T0 >
void operation::extract_template(void (*)( T0& ))
{
    this->profiles[0] = &breakdown(*static_cast<T0*>(this->arguments[0].get()));
}
template < typename T0 >
void operation::prototype_template(void (*)( T0& ))
{
    ((void (*)( T0& ))this->operation_ptr)
    ( *static_cast<T0*>(this->arguments[0].get()) );
}
template < typename T0 > void operation::extract_template(void (*)( pinned T0& )) { this->profiles[0] = &breakdown(*static_cast<T0*>(this->arguments[0].get())); this->pin = this->profiles[0]; } template < typename T0 > void operation::prototype_template(void (*)( pinned T0& )) { ((void (*)( pinned T0& ))this->operation_ptr) ( marked *static_cast<T0*>(this->arguments[0].get()) ); }

template< typename FP, typename T0 , typename T1 >
operation::operation( FP op, T0 *arg0 , T1 *arg1 ){
    this->init();
    this->operation_ptr = (void(*)())op;
    this->count = 2;
    this->profiles = (p_profile**)malloc(sizeof(p_profile*)*this->count);
    this->arguments = new boost::shared_ptr<void>[this->count];
    this->arguments[0] = boost::shared_ptr<void>(get_handle(*arg0)); this->arguments[1] = boost::shared_ptr<void>(get_handle(*arg1));
    void(operation::*ptr)(FP); ptr = &operation::prototype_template;
    this->prototype = (void(operation::*)())ptr;
    ptr = &operation::extract_template;
    this->extract = (void(operation::*)())ptr;
}
template < typename T0 , typename T1 >
void operation::extract_template(void (*)( T0& , T1& ))
{
    this->profiles[0] = &breakdown(*static_cast<T0*>(this->arguments[0].get())); this->profiles[1] = &breakdown(*static_cast<T1*>(this->arguments[1].get()));
}
template < typename T0 , typename T1 >
void operation::prototype_template(void (*)( T0& , T1& ))
{
    ((void (*)( T0& , T1& ))this->operation_ptr)
    ( *static_cast<T0*>(this->arguments[0].get()) , *static_cast<T1*>(this->arguments[1].get()) );
}
template < typename T0 , typename T1 > void operation::extract_template(void (*)( pinned T0& , T1& )) { this->profiles[0] = &breakdown(*static_cast<T0*>(this->arguments[0].get())); this->profiles[1] = &breakdown(*static_cast<T1*>(this->arguments[1].get())); this->pin = this->profiles[0]; } template < typename T0 , typename T1 > void operation::prototype_template(void (*)( pinned T0& , T1& )) { ((void (*)( pinned T0& , T1& ))this->operation_ptr) ( marked *static_cast<T0*>(this->arguments[0].get()) , *static_cast<T1*>(this->arguments[1].get()) ); } template < typename T0 , typename T1 > void operation::extract_template(void (*)( T0& , pinned T1& )) { this->profiles[0] = &breakdown(*static_cast<T0*>(this->arguments[0].get())); this->profiles[1] = &breakdown(*static_cast<T1*>(this->arguments[1].get())); this->pin = this->profiles[1]; } template < typename T0 , typename T1 > void operation::prototype_template(void (*)( T0& , pinned T1& )) { ((void (*)( T0& , pinned T1& ))this->operation_ptr) ( *static_cast<T0*>(this->arguments[0].get()) , marked *static_cast<T1*>(this->arguments[1].get()) ); }

template< typename FP, typename T0 , typename T1 , typename T2 >
operation::operation( FP op, T0 *arg0 , T1 *arg1 , T2 *arg2 ){
    this->init();
    this->operation_ptr = (void(*)())op;
    this->count = 3;
    this->profiles = (p_profile**)malloc(sizeof(p_profile*)*this->count);
    this->arguments = new boost::shared_ptr<void>[this->count];
    this->arguments[0] = boost::shared_ptr<void>(get_handle(*arg0)); this->arguments[1] = boost::shared_ptr<void>(get_handle(*arg1)); this->arguments[2] = boost::shared_ptr<void>(get_handle(*arg2));
    void(operation::*ptr)(FP); ptr = &operation::prototype_template;
    this->prototype = (void(operation::*)())ptr;
    ptr = &operation::extract_template;
    this->extract = (void(operation::*)())ptr;
}
template < typename T0 , typename T1 , typename T2 >
void operation::extract_template(void (*)( T0& , T1& , T2& ))
{
    this->profiles[0] = &breakdown(*static_cast<T0*>(this->arguments[0].get())); this->profiles[1] = &breakdown(*static_cast<T1*>(this->arguments[1].get())); this->profiles[2] = &breakdown(*static_cast<T2*>(this->arguments[2].get()));
}
template < typename T0 , typename T1 , typename T2 >
void operation::prototype_template(void (*)( T0& , T1& , T2& ))
{
    ((void (*)( T0& , T1& , T2& ))this->operation_ptr)
    ( *static_cast<T0*>(this->arguments[0].get()) , *static_cast<T1*>(this->arguments[1].get()) , *static_cast<T2*>(this->arguments[2].get()) );
}
template < typename T0 , typename T1 , typename T2 > void operation::extract_template(void (*)( pinned T0& , T1& , T2& )) { this->profiles[0] = &breakdown(*static_cast<T0*>(this->arguments[0].get())); this->profiles[1] = &breakdown(*static_cast<T1*>(this->arguments[1].get())); this->profiles[2] = &breakdown(*static_cast<T2*>(this->arguments[2].get())); this->pin = this->profiles[0]; } template < typename T0 , typename T1 , typename T2 > void operation::prototype_template(void (*)( pinned T0& , T1& , T2& )) { ((void (*)( pinned T0& , T1& , T2& ))this->operation_ptr) ( marked *static_cast<T0*>(this->arguments[0].get()) , *static_cast<T1*>(this->arguments[1].get()) , *static_cast<T2*>(this->arguments[2].get()) ); } template < typename T0 , typename T1 , typename T2 > void operation::extract_template(void (*)( T0& , pinned T1& , T2& )) { this->profiles[0] = &breakdown(*static_cast<T0*>(this->arguments[0].get())); this->profiles[1] = &breakdown(*static_cast<T1*>(this->arguments[1].get())); this->profiles[2] = &breakdown(*static_cast<T2*>(this->arguments[2].get())); this->pin = this->profiles[1]; } template < typename T0 , typename T1 , typename T2 > void operation::prototype_template(void (*)( T0& , pinned T1& , T2& )) { ((void (*)( T0& , pinned T1& , T2& ))this->operation_ptr) ( *static_cast<T0*>(this->arguments[0].get()) , marked *static_cast<T1*>(this->arguments[1].get()) , *static_cast<T2*>(this->arguments[2].get()) ); } template < typename T0 , typename T1 , typename T2 > void operation::extract_template(void (*)( T0& , T1& , pinned T2& )) { this->profiles[0] = &breakdown(*static_cast<T0*>(this->arguments[0].get())); this->profiles[1] = &breakdown(*static_cast<T1*>(this->arguments[1].get())); this->profiles[2] = &breakdown(*static_cast<T2*>(this->arguments[2].get())); this->pin = this->profiles[2]; } template < typename T0 , typename T1 , typename T2 > void operation::prototype_template(void (*)( T0& , T1& , pinned T2& )) { ((void (*)( T0& , T1& , pinned T2& ))this->operation_ptr) ( *static_cast<T0*>(this->arguments[0].get()) , *static_cast<T1*>(this->arguments[1].get()) , marked *static_cast<T2*>(this->arguments[2].get()) ); }

template< typename FP, typename T0 , typename T1 , typename T2 , typename T3 >
operation::operation( FP op, T0 *arg0 , T1 *arg1 , T2 *arg2 , T3 *arg3 ){
    this->init();
    this->operation_ptr = (void(*)())op;
    this->count = 4;
    this->profiles = (p_profile**)malloc(sizeof(p_profile*)*this->count);
    this->arguments = new boost::shared_ptr<void>[this->count];
    this->arguments[0] = boost::shared_ptr<void>(get_handle(*arg0)); this->arguments[1] = boost::shared_ptr<void>(get_handle(*arg1)); this->arguments[2] = boost::shared_ptr<void>(get_handle(*arg2)); this->arguments[3] = boost::shared_ptr<void>(get_handle(*arg3));
    void(operation::*ptr)(FP); ptr = &operation::prototype_template;
    this->prototype = (void(operation::*)())ptr;
    ptr = &operation::extract_template;
    this->extract = (void(operation::*)())ptr;
}
template < typename T0 , typename T1 , typename T2 , typename T3 >
void operation::extract_template(void (*)( T0& , T1& , T2& , T3& ))
{
    this->profiles[0] = &breakdown(*static_cast<T0*>(this->arguments[0].get())); this->profiles[1] = &breakdown(*static_cast<T1*>(this->arguments[1].get())); this->profiles[2] = &breakdown(*static_cast<T2*>(this->arguments[2].get())); this->profiles[3] = &breakdown(*static_cast<T3*>(this->arguments[3].get()));
}
template < typename T0 , typename T1 , typename T2 , typename T3 >
void operation::prototype_template(void (*)( T0& , T1& , T2& , T3& ))
{
    ((void (*)( T0& , T1& , T2& , T3& ))this->operation_ptr)
    ( *static_cast<T0*>(this->arguments[0].get()) , *static_cast<T1*>(this->arguments[1].get()) , *static_cast<T2*>(this->arguments[2].get()) , *static_cast<T3*>(this->arguments[3].get()) );
}
template < typename T0 , typename T1 , typename T2 , typename T3 > void operation::extract_template(void (*)( pinned T0& , T1& , T2& , T3& )) { this->profiles[0] = &breakdown(*static_cast<T0*>(this->arguments[0].get())); this->profiles[1] = &breakdown(*static_cast<T1*>(this->arguments[1].get())); this->profiles[2] = &breakdown(*static_cast<T2*>(this->arguments[2].get())); this->profiles[3] = &breakdown(*static_cast<T3*>(this->arguments[3].get())); this->pin = this->profiles[0]; } template < typename T0 , typename T1 , typename T2 , typename T3 > void operation::prototype_template(void (*)( pinned T0& , T1& , T2& , T3& )) { ((void (*)( pinned T0& , T1& , T2& , T3& ))this->operation_ptr) ( marked *static_cast<T0*>(this->arguments[0].get()) , *static_cast<T1*>(this->arguments[1].get()) , *static_cast<T2*>(this->arguments[2].get()) , *static_cast<T3*>(this->arguments[3].get()) ); } template < typename T0 , typename T1 , typename T2 , typename T3 > void operation::extract_template(void (*)( T0& , pinned T1& , T2& , T3& )) { this->profiles[0] = &breakdown(*static_cast<T0*>(this->arguments[0].get())); this->profiles[1] = &breakdown(*static_cast<T1*>(this->arguments[1].get())); this->profiles[2] = &breakdown(*static_cast<T2*>(this->arguments[2].get())); this->profiles[3] = &breakdown(*static_cast<T3*>(this->arguments[3].get())); this->pin = this->profiles[1]; } template < typename T0 , typename T1 , typename T2 , typename T3 > void operation::prototype_template(void (*)( T0& , pinned T1& , T2& , T3& )) { ((void (*)( T0& , pinned T1& , T2& , T3& ))this->operation_ptr) ( *static_cast<T0*>(this->arguments[0].get()) , marked *static_cast<T1*>(this->arguments[1].get()) , *static_cast<T2*>(this->arguments[2].get()) , *static_cast<T3*>(this->arguments[3].get()) ); } template < typename T0 , typename T1 , typename T2 , typename T3 > void operation::extract_template(void (*)( T0& , T1& , pinned T2& , T3& )) { this->profiles[0] = &breakdown(*static_cast<T0*>(this->arguments[0].get())); this->profiles[1] = &breakdown(*static_cast<T1*>(this->arguments[1].get())); this->profiles[2] = &breakdown(*static_cast<T2*>(this->arguments[2].get())); this->profiles[3] = &breakdown(*static_cast<T3*>(this->arguments[3].get())); this->pin = this->profiles[2]; } template < typename T0 , typename T1 , typename T2 , typename T3 > void operation::prototype_template(void (*)( T0& , T1& , pinned T2& , T3& )) { ((void (*)( T0& , T1& , pinned T2& , T3& ))this->operation_ptr) ( *static_cast<T0*>(this->arguments[0].get()) , *static_cast<T1*>(this->arguments[1].get()) , marked *static_cast<T2*>(this->arguments[2].get()) , *static_cast<T3*>(this->arguments[3].get()) ); } template < typename T0 , typename T1 , typename T2 , typename T3 > void operation::extract_template(void (*)( T0& , T1& , T2& , pinned T3& )) { this->profiles[0] = &breakdown(*static_cast<T0*>(this->arguments[0].get())); this->profiles[1] = &breakdown(*static_cast<T1*>(this->arguments[1].get())); this->profiles[2] = &breakdown(*static_cast<T2*>(this->arguments[2].get())); this->profiles[3] = &breakdown(*static_cast<T3*>(this->arguments[3].get())); this->pin = this->profiles[3]; } template < typename T0 , typename T1 , typename T2 , typename T3 > void operation::prototype_template(void (*)( T0& , T1& , T2& , pinned T3& )) { ((void (*)( T0& , T1& , T2& , pinned T3& ))this->operation_ptr) ( *static_cast<T0*>(this->arguments[0].get()) , *static_cast<T1*>(this->arguments[1].get()) , *static_cast<T2*>(this->arguments[2].get()) , marked *static_cast<T3*>(this->arguments[3].get()) ); }

template< typename FP, typename T0 , typename T1 , typename T2 , typename T3 , typename T4 >
operation::operation( FP op, T0 *arg0 , T1 *arg1 , T2 *arg2 , T3 *arg3 , T4 *arg4 ){
    this->init();
    this->operation_ptr = (void(*)())op;
    this->count = 5;
    this->profiles = (p_profile**)malloc(sizeof(p_profile*)*this->count);
    this->arguments = new boost::shared_ptr<void>[this->count];
    this->arguments[0] = boost::shared_ptr<void>(get_handle(*arg0)); this->arguments[1] = boost::shared_ptr<void>(get_handle(*arg1)); this->arguments[2] = boost::shared_ptr<void>(get_handle(*arg2)); this->arguments[3] = boost::shared_ptr<void>(get_handle(*arg3)); this->arguments[4] = boost::shared_ptr<void>(get_handle(*arg4));
    void(operation::*ptr)(FP); ptr = &operation::prototype_template;
    this->prototype = (void(operation::*)())ptr;
    ptr = &operation::extract_template;
    this->extract = (void(operation::*)())ptr;
}
template < typename T0 , typename T1 , typename T2 , typename T3 , typename T4 >
void operation::extract_template(void (*)( T0& , T1& , T2& , T3& , T4& ))
{
    this->profiles[0] = &breakdown(*static_cast<T0*>(this->arguments[0].get())); this->profiles[1] = &breakdown(*static_cast<T1*>(this->arguments[1].get())); this->profiles[2] = &breakdown(*static_cast<T2*>(this->arguments[2].get())); this->profiles[3] = &breakdown(*static_cast<T3*>(this->arguments[3].get())); this->profiles[4] = &breakdown(*static_cast<T4*>(this->arguments[4].get()));
}
template < typename T0 , typename T1 , typename T2 , typename T3 , typename T4 >
void operation::prototype_template(void (*)( T0& , T1& , T2& , T3& , T4& ))
{
    ((void (*)( T0& , T1& , T2& , T3& , T4& ))this->operation_ptr)
    ( *static_cast<T0*>(this->arguments[0].get()) , *static_cast<T1*>(this->arguments[1].get()) , *static_cast<T2*>(this->arguments[2].get()) , *static_cast<T3*>(this->arguments[3].get()) , *static_cast<T4*>(this->arguments[4].get()) );
}
template < typename T0 , typename T1 , typename T2 , typename T3 , typename T4 > void operation::extract_template(void (*)( pinned T0& , T1& , T2& , T3& , T4& )) { this->profiles[0] = &breakdown(*static_cast<T0*>(this->arguments[0].get())); this->profiles[1] = &breakdown(*static_cast<T1*>(this->arguments[1].get())); this->profiles[2] = &breakdown(*static_cast<T2*>(this->arguments[2].get())); this->profiles[3] = &breakdown(*static_cast<T3*>(this->arguments[3].get())); this->profiles[4] = &breakdown(*static_cast<T4*>(this->arguments[4].get())); this->pin = this->profiles[0]; } template < typename T0 , typename T1 , typename T2 , typename T3 , typename T4 > void operation::prototype_template(void (*)( pinned T0& , T1& , T2& , T3& , T4& )) { ((void (*)( pinned T0& , T1& , T2& , T3& , T4& ))this->operation_ptr) ( marked *static_cast<T0*>(this->arguments[0].get()) , *static_cast<T1*>(this->arguments[1].get()) , *static_cast<T2*>(this->arguments[2].get()) , *static_cast<T3*>(this->arguments[3].get()) , *static_cast<T4*>(this->arguments[4].get()) ); } template < typename T0 , typename T1 , typename T2 , typename T3 , typename T4 > void operation::extract_template(void (*)( T0& , pinned T1& , T2& , T3& , T4& )) { this->profiles[0] = &breakdown(*static_cast<T0*>(this->arguments[0].get())); this->profiles[1] = &breakdown(*static_cast<T1*>(this->arguments[1].get())); this->profiles[2] = &breakdown(*static_cast<T2*>(this->arguments[2].get())); this->profiles[3] = &breakdown(*static_cast<T3*>(this->arguments[3].get())); this->profiles[4] = &breakdown(*static_cast<T4*>(this->arguments[4].get())); this->pin = this->profiles[1]; } template < typename T0 , typename T1 , typename T2 , typename T3 , typename T4 > void operation::prototype_template(void (*)( T0& , pinned T1& , T2& , T3& , T4& )) { ((void (*)( T0& , pinned T1& , T2& , T3& , T4& ))this->operation_ptr) ( *static_cast<T0*>(this->arguments[0].get()) , marked *static_cast<T1*>(this->arguments[1].get()) , *static_cast<T2*>(this->arguments[2].get()) , *static_cast<T3*>(this->arguments[3].get()) , *static_cast<T4*>(this->arguments[4].get()) ); } template < typename T0 , typename T1 , typename T2 , typename T3 , typename T4 > void operation::extract_template(void (*)( T0& , T1& , pinned T2& , T3& , T4& )) { this->profiles[0] = &breakdown(*static_cast<T0*>(this->arguments[0].get())); this->profiles[1] = &breakdown(*static_cast<T1*>(this->arguments[1].get())); this->profiles[2] = &breakdown(*static_cast<T2*>(this->arguments[2].get())); this->profiles[3] = &breakdown(*static_cast<T3*>(this->arguments[3].get())); this->profiles[4] = &breakdown(*static_cast<T4*>(this->arguments[4].get())); this->pin = this->profiles[2]; } template < typename T0 , typename T1 , typename T2 , typename T3 , typename T4 > void operation::prototype_template(void (*)( T0& , T1& , pinned T2& , T3& , T4& )) { ((void (*)( T0& , T1& , pinned T2& , T3& , T4& ))this->operation_ptr) ( *static_cast<T0*>(this->arguments[0].get()) , *static_cast<T1*>(this->arguments[1].get()) , marked *static_cast<T2*>(this->arguments[2].get()) , *static_cast<T3*>(this->arguments[3].get()) , *static_cast<T4*>(this->arguments[4].get()) ); } template < typename T0 , typename T1 , typename T2 , typename T3 , typename T4 > void operation::extract_template(void (*)( T0& , T1& , T2& , pinned T3& , T4& )) { this->profiles[0] = &breakdown(*static_cast<T0*>(this->arguments[0].get())); this->profiles[1] = &breakdown(*static_cast<T1*>(this->arguments[1].get())); this->profiles[2] = &breakdown(*static_cast<T2*>(this->arguments[2].get())); this->profiles[3] = &breakdown(*static_cast<T3*>(this->arguments[3].get())); this->profiles[4] = &breakdown(*static_cast<T4*>(this->arguments[4].get())); this->pin = this->profiles[3]; } template < typename T0 , typename T1 , typename T2 , typename T3 , typename T4 > void operation::prototype_template(void (*)( T0& , T1& , T2& , pinned T3& , T4& )) { ((void (*)( T0& , T1& , T2& , pinned T3& , T4& ))this->operation_ptr) ( *static_cast<T0*>(this->arguments[0].get()) , *static_cast<T1*>(this->arguments[1].get()) , *static_cast<T2*>(this->arguments[2].get()) , marked *static_cast<T3*>(this->arguments[3].get()) , *static_cast<T4*>(this->arguments[4].get()) ); } template < typename T0 , typename T1 , typename T2 , typename T3 , typename T4 > void operation::extract_template(void (*)( T0& , T1& , T2& , T3& , pinned T4& )) { this->profiles[0] = &breakdown(*static_cast<T0*>(this->arguments[0].get())); this->profiles[1] = &breakdown(*static_cast<T1*>(this->arguments[1].get())); this->profiles[2] = &breakdown(*static_cast<T2*>(this->arguments[2].get())); this->profiles[3] = &breakdown(*static_cast<T3*>(this->arguments[3].get())); this->profiles[4] = &breakdown(*static_cast<T4*>(this->arguments[4].get())); this->pin = this->profiles[4]; } template < typename T0 , typename T1 , typename T2 , typename T3 , typename T4 > void operation::prototype_template(void (*)( T0& , T1& , T2& , T3& , pinned T4& )) { ((void (*)( T0& , T1& , T2& , T3& , pinned T4& ))this->operation_ptr) ( *static_cast<T0*>(this->arguments[0].get()) , *static_cast<T1*>(this->arguments[1].get()) , *static_cast<T2*>(this->arguments[2].get()) , *static_cast<T3*>(this->arguments[3].get()) , marked *static_cast<T4*>(this->arguments[4].get()) ); }
}
