// nested inside ambient.hpp in ambient namespace

void memoryfence(){
    ambient::spin();
    //scope.get_group()->get_manager()->spin_loop(); // just spin for now
}

template <typename FL, typename FC, class T0, class T1, class T2>
void push(FL l_kernel, FC c_kernel, T0& arg0, T1& arg1, T2& arg2){
    if(breakdown(arg0).state == PROXY) pin(arg0, arg2);
    if(breakdown(arg1).state == PROXY) pin(arg1, arg2);
    ambient::engine.push(new core::operation(l_kernel, &arg0, &arg1, &arg2),
			 new core::operation(c_kernel, &arg0, &arg1, &arg2));
}
template <typename ST, typename FL, typename FC, class T0, class T1>
ST& push(FL l_kernel, FC c_kernel, T0& arg0, T1& arg1){
    void_pt* handle = new void_pt((ST*)NULL);
    ST out(handle);
    push(l_kernel, c_kernel, arg0, arg1, out);
    return *(ST*)out.self; // trying to omit casting copy
}
#include "ambient/interface/push.pp.hpp" // all other variants of push

template <typename L, typename R>
void pin(L& proxy_object, const R& real_object){
    void_pt& proxy = breakdown(proxy_object);
    void_pt& real  = breakdown(real_object);
    ((livelong<R,REPLICA>*)(get_handle(real_object).get()))->use_count += 2; // avoiding false deallocation (which is memory leak)
    proxy.profile   = real.profile;
    proxy.set_dim(real.get_dim());
    real.imitate(&proxy); // copy proxy settings to the actual profile
}
template<typename T>
void assign(const T& ref, int i, int j)
{
    void_pt& profile = breakdown(ref);
    profile.layout->request(i, j);
}
template<typename T>
void assign(T& ref, int i, int j)
{
    void_pt& profile = breakdown(ref);
    profile.layout->record(i, j);
}
template<typename T>
inline std::pair<unsigned int*,size_t> get_id(T& ref)
{
    return breakdown(ref).get_id();
}
template<typename T>
inline dim2 get_dim(T& ref)
{
    return breakdown(ref).get_dim();
}
template<typename T>
inline dim2 get_distr_dim(T& ref)
{
    return breakdown(ref).get_distr_dim();
}
template<typename T>
inline dim2 get_gpu_dim(T& ref)
{
    return breakdown(ref).get_gpu_dim();
}
template<typename T>
inline dim2 get_grid_dim(T& ref)
{
    return breakdown(ref).get_grid_dim();
}
template<typename T>
inline dim2 get_group_dim(T& ref)
{
    return breakdown(ref).get_group_dim();
}
template<typename T>
inline dim2 get_group_t_dim(T& ref)
{
    return breakdown(ref).get_group_t_dim();
}
template<typename T>
inline dim2 get_item_dim(T& ref)
{
    return breakdown(ref).get_item_dim();
}
template<typename T>
inline dim2 get_group_id(T& ref)
{
    return breakdown(ref).get_group_id();
}

