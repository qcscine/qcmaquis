template<typename T> 
void info(T& obj){
    if(rank.is_master(scope.get_group())){
        void_pt& p = breakdown(obj);
        printf("Matrix %d:%d size of the task is %d x %d groups sized %d x %d items of %d x %d elements\n", 
               *p.group_id, p.id, p.get_grid_dim().y, p.get_grid_dim().x, p.get_group_dim().y, p.get_group_dim().x, p.get_item_dim().x, p.get_item_dim().y);
    }
}

template<typename T>
void block_2d_cycle_assign(T& target)
{
///////////////////////////////////////////// 2D-block-cyclic decomposition
    int np = 1; // can be a function arg   // process grid's num of rows 
    int nq = (int)(scope.get_size() / np); // process grid's num of cols 
    int rank_i = (int)(scope.get_rank() / nq); // process row
    int rank_j = (int)(scope.get_rank() % nq); // process col
///////////////////////////////////////////////////////////////////////////
    for(int i = rank_i; i < get_grid_dim(target).y; i += np){
        for(int j = rank_j; j < get_grid_dim(target).x; j += nq){
            assign(target, i, j);
        }
    }
}

void gemm_l_kernel(pinned const p_dense_matrix<double>& a, const p_dense_matrix<double>& b, p_dense_matrix<double>& c)
{
    breakdown(a) >> dim3(), dim3(), dim3();
    scope_select("2 from ambient as gemm where master is 0"); // todo: correct the naming issue
    if(!scope.involved()) return;

    zout << "2d-block-cyclic decomposition kernel in gemm ("<< ambient::rank() <<"):\n"; info(a); info(b); info(c);

    block_2d_cycle_assign(a);
    block_2d_cycle_assign(b);
    block_2d_cycle_assign(c);
}

void copy_l_kernel(p_dense_matrix<double>& a_copy, pinned const p_dense_matrix<double>& a)
{
    scope_select("2 from ambient as copy_ground where master is 0");
    if(!scope.involved()) return; // out of scope quick exit
    zout << "2d-block-cyclic decomposition kernel in copy ("<< ambient::rank() <<"):\n"; info(a); info(a_copy);

    block_2d_cycle_assign(a);
    block_2d_cycle_assign(a_copy);
}

void mem_bound_l_kernel(const p_dense_matrix<double>& a, const p_dense_matrix<double>& b, pinned p_dense_matrix<double>& c)
{
    scope_select("2 from ambient as work where master is 0");
    scope_retain("2 from ambient as work_storage");
    if(!scope.involved()) return; // out of scope quick exit
    zout << "2d-block-cyclic decomposition kernel in membound ("<< ambient::rank() <<"):\n"; info(a); info(b); info(c);

    block_2d_cycle_assign(a);
    block_2d_cycle_assign(b);
    block_2d_cycle_assign(c);
}

void scale_l_kernel(const p_dense_matrix<double>& m, const double& t, pinned p_dense_matrix<double>& out){}

/////////////////////
// testing kernels // 

void single_integer_l_kernel(int& input){
    scope_select("* from ambient as single_integer_work where master is 0");
    if(!scope.involved()) return;
    zout << "single integer kernel: input is " << input << "\n";
    input = 0;
}
