// nested inside ambient.hpp in ambient namespace

template<class T>
void matrix_i_kernel(memblock* grp)
{
    int n = grp->get_group_t_dim().x;
    int m = grp->get_group_t_dim().y;
    int ld = m;

    //memset(grp->data,0,m*n*sizeof(T));
    for(int j=0; j<n; j++)
    for(int i=0; i<m; i++){
        ((double*)grp->data)[j*ld+i] = drand48();
    }
}

template<class T>
void nullify(memblock* grp)
{
    int n = grp->get_group_t_dim().x;
    int m = grp->get_group_t_dim().y;

    memset(grp->data,0,m*n*sizeof(T));
}
