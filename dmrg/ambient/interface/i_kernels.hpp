// nested inside ambient.hpp in ambient namespace
using namespace blas;

void matrix_i_kernel(workgroup* grp){
    // real temporary initialization
    int n = grp->get_group_dim().x*grp->get_item_dim().x;
    int m = grp->get_group_dim().y*grp->get_item_dim().y;
    int ld = m;

    for(int j=0; j<n; j++)
    for(int i=0; i<m; i++){
        ((double*)grp->data)[j*ld+i] = i+j;
    }
}
