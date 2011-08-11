// nested inside ambient.hpp in ambient namespace

void memoryfence(){
    ambient::spin();
    //scope.get_group()->spin_loop(); // just spin for now
}

#include "ambient/interface/push.pp.hpp" // all variants of push

template<typename T>
void assign(const T& ref, int i, int j)
{
    void_pt& profile = breakdown(ref);
    dim2 work_blocks(profile.get_work_dim().x / profile.get_mem_dim().x,
                     profile.get_work_dim().y / profile.get_mem_dim().y);
    int ii = i*work_blocks.y;
    int jj = j*work_blocks.x;

// be careful not to overflow the borders!
    for(int i = ii; i < ii+work_blocks.y; i++)
        for(int j = jj; j < jj+work_blocks.x; j++)
            profile.layout->request(i, j);
}
template<typename T>
void assign(T& ref, int i, int j)
{
    void_pt& profile = breakdown(ref);
    dim2 work_blocks(profile.get_work_dim().x / profile.get_mem_dim().x,
                     profile.get_work_dim().y / profile.get_mem_dim().y);
    int ii = i*work_blocks.y;
    int jj = j*work_blocks.x;

// be careful not to overflow the borders!
    for(int i = ii; i < ii+work_blocks.y; i++)
        for(int j = jj; j < jj+work_blocks.x; j++)
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
inline dim2 get_work_dim(T& ref)
{
    return breakdown(ref).get_work_dim();
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
inline dim2 get_mem_dim(T& ref)
{
    return breakdown(ref).get_mem_dim();
}
template<typename T>
inline dim2 get_mem_t_dim(T& ref)
{
    return breakdown(ref).get_mem_t_dim();
}
template<typename T>
inline dim2 get_item_dim(T& ref)
{
    return breakdown(ref).get_item_dim();
}
template<typename T>
inline dim2 get_block_id(T& ref)
{
    return breakdown(ref).get_block_id();
}

