template<typename T> 
void info(T& obj){
    if(rank.is_master(scope.get_group())){
        void_pt& p = breakdown(obj);
        printf("Matrix %d:%d size of the task is %d x %d groups sized %d x %d items of %d x %d elements\n", 
               *p.group_id, p.id, p.get_grid_dim().y, p.get_grid_dim().x, p.get_group_dim().y, p.get_group_dim().x, p.get_item_dim().x, p.get_item_dim().y);
    }
}

void mem_bound_l_kernel(const p_dense_matrix<double>& a, const p_dense_matrix<double>& b, pinned p_dense_matrix<double>& out){
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

void gemm_l_kernel(const p_dense_matrix<double>& a, const p_dense_matrix<double>& b, pinned p_dense_matrix<double>& out)
{
    scope_select("0.5 from ambient as work where master is 0"); // todo: correct the naming issue
    if(!scope.involved()) return;

    zout << "2d-block-cyclic decomposition kernel in gemm:\n"; info(a); info(b); info(out);
///////////////////////////////////////////// 2D-block-cyclic decomposition
    int np = 1; // can be a function arg   // process grid's num of rows 
    int nq = (int)(scope.get_size() / np); // process grid's num of cols 
    int rank_i = (int)(scope.get_rank() / nq); // process row
    int rank_j = (int)(scope.get_rank() % nq); // process col
///////////////////////////////////////////////////////////////////////////
    for(int i = rank_i; i < get_grid_dim(out).y; i += np){
        for(int j = rank_j; j < get_grid_dim(out).x; j += nq){
            assign(a,   i, 0);
            assign(b,   0, j);
            assign(out, i, j);
        }
    }
}

void scale_l_kernel(const p_dense_matrix<double>& m, const double& t, pinned p_dense_matrix<double>& out){
// todo
}

void pdgemm_l_kernel(p_dense_matrix<double>& a, p_dense_matrix<double>& b, p_dense_matrix<double>& c){
    scope_select("* from ambient as work_redist where master is 0"); // fix the "same name" problem!
    if(!scope.involved()) return;

    zout << "2d-block-cyclic decomposition kernel:\n"; info(a); info(b); info(c);
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
    for(int i = rank_i; i < get_grid_dim(a).y; i += np){
        for(int j = rank_j; j < get_grid_dim(a).x; j += nq){
            assign(b, i, j);
        }
    }
    for(int i = rank_i; i < get_grid_dim(a).y; i += np){
        for(int j = rank_j; j < get_grid_dim(a).x; j += nq){
            assign(c, i, j);
        }
    }
}

void single_integer_l_kernel(int& input){
    scope_select("* from ambient as single_integer_work where master is 0");
    if(!scope.involved()) return;
    zout << "single integer kernel: input is " << input << "\n";
    input = 0;
}

