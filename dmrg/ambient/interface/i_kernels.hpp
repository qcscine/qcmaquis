// nested inside ambient.hpp in ambient namespace
using namespace blas;

void matrix_i_kernel(workgroup* grp){
    // dumb 0-initialization
    memset(grp->data, 0, grp->get_group_dim().y*grp->get_item_dim().y*grp->get_group_dim().x*grp->get_item_dim().x*grp->get_profile()->type_size);
    // debug fillidfdsfds:
    size_t size = grp->get_group_dim().y*grp->get_item_dim().y*grp->get_group_dim().x*grp->get_item_dim().x;
    for(size_t i=0; i < size; i++)
        ((double*)grp->data)[i] = i;
}
