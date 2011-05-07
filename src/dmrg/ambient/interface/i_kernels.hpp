// nested inside ambient.hpp in ambient namespace

template<typename T>
void random_i_kernel(pinned p_dense_matrix<T>& a)
{
    size_t i = get_block_id(a).y;
    size_t j = get_block_id(a).x;
    size_t m = get_mem_t_dim(a).y;
    size_t n = get_mem_t_dim(a).x;
    size_t ld = m;
    size_t sd = n;
    T* ad = current(a)(i,j);

    if(i == (get_grid_dim(a).y-1) && get_dim(a).y%get_mem_t_dim(a).y != 0) 
        m = get_dim(a).y % get_mem_t_dim(a).y;
    if(j == (get_grid_dim(a).x-1) && get_dim(a).x%get_mem_t_dim(a).x != 0) 
        n = get_dim(a).x % get_mem_t_dim(a).x;
    
    if(m != ld || n != sd){
        for(size_t jj=0; jj<n; jj++){
            for(size_t ii=0; ii<m; ii++) ad[jj*ld+ii] = drand48();
            if(m != ld) memset(&ad[jj*ld+m], 0, sizeof(T)*(ld-m));
        }
        if(n != sd) 
            memset(&ad[n*ld], 0, sizeof(T)*(sd-n)*ld);
    }else
        for(size_t ii=0; ii<n*m; ii++) ad[ii] = drand48();
}

template<typename T>
void nullcut_i_kernel(pinned p_dense_matrix<T>& a)
{
    size_t i = get_block_id(a).y;
    size_t j = get_block_id(a).x;
    size_t m = get_mem_t_dim(a).y;
    size_t n = get_mem_t_dim(a).x;
    T* ad = current(a)(i,j);
    memset(ad, 0, m*n*sizeof(T));
}

template<class T>
void matrix_i_kernel(memblock* grp)
{
    int n = grp->get_mem_t_dim().x;
    int m = grp->get_mem_t_dim().y;
    int ld = m;

    //memset(grp->data,0,m*n*sizeof(T));
    for(int j=0; j<n; j++)
    for(int i=0; i<m; i++){
        ((T*)grp->data)[j*ld+i] = drand48();
    }
}


template<class T>
void nullify(memblock* grp)
{
    int n = grp->get_mem_t_dim().x;
    int m = grp->get_mem_t_dim().y;

    memset(grp->data,0,m*n*sizeof(T));
}
