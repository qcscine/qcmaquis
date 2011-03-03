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

template<typename T> 
void info(T& obj){
    if(rank.is_master(scope.get_group())){
        void_pt& p = breakdown(obj);
        printf("Matrix %d:%d size of the task is %d x %d groups sized %d x %d items of %d x %d elements\n", 
               *p.group_id, p.id, p.get_grid_dim().y, p.get_grid_dim().x, p.get_group_dim().y, p.get_group_dim().x, p.get_item_dim().x, p.get_item_dim().y);
    }
}

void plus_l_kernel(const p_dense_matrix<double>& a, const p_dense_matrix<double>& b, pinned p_dense_matrix<double>& out){
//    a >> dim3(10,5), dim3(1,1), dim3(10,1); <- kinda non-trivial - need to think
    scope_select("2 from ambient as work where master is 0");
    scope_retain("2 from ambient as work_storage");
    if(!scope.involved()) return; // out of scope quick exit

    info(a); info(b); info(out);

    for(int i=0; i < get_grid_dim(out).y; i++)
        for(int j=0; j < get_grid_dim(out).x; j++)
            if(j % scope.get_size() == scope.get_rank()){
                assign(a,   i, j);
                assign(b,   i, j);
                assign(out, i, j);
            }
}

void plus_c_kernel(const p_dense_matrix<double>& a, const p_dense_matrix<double>& b, pinned p_dense_matrix<double>& out){
    void_pt& profile = breakdown(out);
    double* ad = breakdown(a)(get_group_id(out).x, get_group_id(out).y);
    double* bd = breakdown(b)(get_group_id(out).x, get_group_id(out).y);
    int size = get_group_dim(out).x*get_item_dim(out).x*
               get_group_dim(out).y*get_item_dim(out).y;
//    printf("R%d: Executing plus computation kernel (%d ops)... for out grp %d %d\n", scope.get_rank(), size, get_group_id(out).x, get_group_id(out).y);
//    for(int i=0; i < size; i++){
//        output[i] = ad[i]+bd[i];
//    }
}

void redistribution_l_kernel(p_dense_matrix<double>& a){
    scope_select("* from ambient as work_redist where master is 0");
    if(!scope.involved()) return;

    zout << "2d-block-cyclic decomposition kernel:\n"; info(a);
///////////////////////////////////////////// 2D-block-cyclic decomposition
    int np = 1; // can be a function arg   // process grid's num of rows 
    int nq = (int)(scope.get_size() / np); // process grid's num of cols 
    int rank_i = (int)(scope.get_rank() / nq); // process row
    int rank_j = (int)(scope.get_rank() % nq); // process col
///////////////////////////////////////////////////////////////////////////
    for(int i = rank_i; i < get_grid_dim(a).y; i += np){
        for(int j = rank_j; j < get_grid_dim(a).x; j += nq){
            assign(a, i, j);
        }
    }
}

void redistribution_c_kernel(p_dense_matrix<double>& a){

}

void single_integer_l_kernel(int& input){
    scope_select("* from ambient as single_integer_work where master is 0");
    if(!scope.involved()) return;
    zout << "single integer kernel: input is " << input << "\n";
    input = 0;
}

void single_integer_c_kernel(int& input){
    input += 13;
    zout << "single integer kernel: output is " << input << "\n";
}

