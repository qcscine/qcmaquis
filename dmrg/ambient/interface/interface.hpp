// nested inside ambient.hpp in ambient namespace

template <typename FL, typename FC, class T0, class T1, class T2>
void push(FL l_kernel, FC c_kernel, T0& arg0, T1& arg1, T2& arg2){
    if(breakdown(arg0).proxy) pin(arg0, arg2);
    if(breakdown(arg1).proxy) pin(arg1, arg2);
    ambient::engine.push(new core::operation(l_kernel, &arg0, &arg1, &arg2),
			 new core::operation(c_kernel, &arg0, &arg1, &arg2));
}
template <typename ST, typename FL, typename FC, class T0, class T1>
ST push(FL l_kernel, FC c_kernel, T0& arg0, T1& arg1){
    void_pt* handle = new void_pt((ST*)NULL);
    ST out(handle);
    push(l_kernel, c_kernel, arg0, arg1, out);
    return out;
}
#include "ambient/interface/push.pp.hpp" // all other variants of push

template <typename L, typename R>
void pin(L& proxy_object, const R& real_object){
    void_pt& proxy = breakdown(proxy_object);
    void_pt& real  = breakdown(real_object);
    proxy.profile   = real.profile;
    proxy.set_dim(real.get_dim());
    real.imitate(&proxy); // copy proxy settings to the actual profile
}
template<typename T>
void assign(T& ref, int i, int j, int k)
{
// need to check uniqness here...
    void_pt& profile = breakdown(ref);
    workgroup* group = profile.group(i,j,k);
        printf("%s: p%d: I've accepted group %d %d of id%d\n", scope.get_group()->name, scope.get_rank(), group->i, group->j, (*(group->profile))->id );
    group->owner = ambient::rank();
    profile.layout->update_map_entry(ambient::rank(), i, j, k); // or add_segment_entry
}
template<typename T>
void assign(const T& ref, int i, int j, int k)
{
    void_pt& profile = breakdown(ref);
    workgroup* group = profile.group(i,j,k);
//        printf("%s: p%d: I've accepted group %d %d of id%d\n", scope.get_group()->name, scope.get_rank(), group->i, group->j, (*(group->profile))->id );
    group->owner = ambient::rank();
    profile.layout->update_map_entry(ambient::rank(), i, j, k); // or add_segment_entry
}
template<typename T>
dim3 get_dim(T& ref)
{
    return breakdown(ref).get_dim();
}
template<typename T>
dim3 get_distr_dim(T& ref)
{
    return breakdown(ref).get_distr_dim();
}
template<typename T>
dim3 get_gpu_dim(T& ref)
{
    return breakdown(ref).get_gpu_dim();
}
template<typename T>
dim3 get_grid_dim(T& ref)
{
    return breakdown(ref).get_grid_dim();
}
template<typename T>
dim3 get_group_dim(T& ref)
{
    return breakdown(ref).get_group_dim();
}
template<typename T>
dim3 get_item_dim(T& ref)
{
    return breakdown(ref).get_item_dim();
}
template<typename T>
dim3 get_group_id(T& ref)
{
    return breakdown(ref).get_group_id();
}

