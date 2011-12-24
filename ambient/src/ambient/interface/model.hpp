// nested inside ambient.hpp in ambient namespace


template <typename T, policy P>  
void breakdown_model(parallel_t* profile, const p_dense_matrix<T,P>* ptr)
{
    profile->t_size = sizeof(T);
    profile->dim.x  = ptr->num_cols();
    profile->dim.y  = ptr->num_rows();
}

template <typename T, policy P>
void bind_model(const p_dense_matrix<T,P>* ptr)
{
    resize_bind_model(ptr);
} 

template <typename T, policy P>
void copy_bind_model(const p_dense_matrix<T,P>* ptr)
{
    resize_bind_model(ptr);
} 

template <typename T, policy P>
void resize_bind_model(const p_dense_matrix<T,P>* ptr)
{
    ptr->breakdown()->set_dim(dim2((unsigned int)ptr->num_cols(), (unsigned int)ptr->num_rows()));
}

template <typename T, policy P>
void resize_bind_model(p_dense_matrix<T,P>* ptr, size_t rows, size_t cols)
{
    ptr->latch_meta(); // for unloose objects only!
    ptr->resize(rows, cols);
}

template <typename T>
parallel_t& breakdown(T& obj){
    static parallel_t* null_profile = new parallel_t(); 
    return *null_profile; 
}

template<>
void plus_reduce< p_dense_matrix<ambient::traits::value_type> >(memblock* grp, void* update){
    ambient::traits::value_type* a = (ambient::traits::value_type*)grp->data;
    ambient::traits::value_type* u = (ambient::traits::value_type*)update;
    int n = grp->get_mem_t_dim().x;
    int m = grp->get_mem_t_dim().y;
    int ld = m;
    std::size_t num = m*n; // C - To vectorize
    for(std::size_t i=0; i<num; ++i) a[i] += u[i];
}

