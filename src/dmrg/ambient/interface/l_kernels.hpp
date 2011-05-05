// AMBIENT USAGE TIPS //
//
// Arguments locality:
//
// when executing the kernel be sure to remember that arguments
// that you are feeding to it are bounded by the life time of 
// C++ scope { }. That is either you need to be sure that the
// kernel will be executed inside current scope or allocate
// them inside heap so that they won't be deleted.
// For example for kernel that takes argument "int*& array"
// one can allocate the argument array the following way:
//
// int*const *array = new int*((int*)malloc(sizeof(int)*LEN));
// for(int i = 0; i < LEN; i++) (*array)[i] = ...;
// ambient::push(... *array);
// 
// Here I'm purposedly underlying that array is a pointer to
// constant pointer to ints - that is I don't want it to be
// deleted. Upon the exit from C++ scope the array ptr will
// be deleted but we would be safe as its value is still valid
// because the int*const was alocated inside heap.

template<typename T> 
void info(T& obj){
    if(rank.is_master(scope.get_group())){
        void_pt& p = breakdown(obj);
        printf("M%d:%d / %d x %d groups (each %d x %d items of %d x %d els)\n", 
               *p.group_id, p.id, p.get_grid_dim().y, p.get_grid_dim().x, p.get_mem_dim().y, p.get_mem_dim().x, p.get_item_dim().x, p.get_item_dim().y);
    }
}

template<typename T>
void block_2d_cycle_assign(T& target)
{
///////////////////////////////////////////// 2D-block-cyclic decomposition
    int np = scope.np = 1; // can be a function arg   // process grid's num of rows 
    int nq = scope.nq = (int)(scope.get_size() / np); // process grid's num of cols 
    int rank_i = (int)(scope.get_rank() / nq); // process row
    int rank_j = (int)(scope.get_rank() % nq); // process col
///////////////////////////////////////////////////////////////////////////
    for(int i = rank_i; i < get_grid_dim(target).y; i += np){
        for(int j = rank_j; j < get_grid_dim(target).x; j += nq){
            assign(target, i, j);
        }
    }
}

template<typename T>
void block_outright_assign(T& target)
{
// this no "distribution" is only needed in side scalapack solver (SVD, ...) where we need a contiguous array (no splitting between proc) for the output (schmidt values ...)
///////////////////////////////////////////////////////////////////////////
    for(int i = 0; i < get_grid_dim(target).y; i++){
        for(int j = 0; j < get_grid_dim(target).x; j++){
            assign(target, i, j);
        }
    }
}


void init_double_l_kernel(pinned p_dense_matrix<double> & a)
{
    int num = 1;//get_grid_dim(a).y;
    scope_select(num+" from ambient as init_double where master is 0"); // todo: correct the naming issue
    if(!scope.involved()) return;
//    zout << "2dbcd for "<< num <<" procs in init_double ("<< ambient::rank() <<"):\n"; info(a);

    block_2d_cycle_assign(a);
}

void gemm_l_kernel(pinned const p_dense_matrix<double>& a, const p_dense_matrix<double>& b, p_dense_matrix<double>& c)
{
    int num = 1;//get_grid_dim(a).y;
    scope_select(num+" from ambient as gemm where master is 0 and breakdown contains "+get_id(a));
    if(!scope.involved()) return;
//    zout << "2dbcd for "<< num <<" procs in gemm ("<< ambient::rank() <<"):\n"; info(a); info(b); info(c);

    block_2d_cycle_assign(a);
    block_2d_cycle_assign(b);
    block_2d_cycle_assign(c);
}

void copy_l_kernel(p_dense_matrix<double>& ac, pinned const p_dense_matrix<double>& a)
{
    scope_select("1 from ambient as copy where master is 0");
    if(!scope.involved()) return;
//    zout << "2dbcd in copy ("<< ambient::rank() <<"):\n"; info(ac); info(a);

    block_2d_cycle_assign(ac);
    block_2d_cycle_assign(a);
}

void resize_l_kernel(p_dense_matrix<double>& a, const size_t& rows, const size_t& cols)
{
    breakdown(a).set_dim(ambient::dim2(cols,rows));
    scope_select("1 from ambient as resize where master is 0");
    if(!scope.involved()) return;
//    zout << "2dbcd in resize ("<< ambient::rank() <<"):\n"; info(a);

    block_2d_cycle_assign(a);
}

void remove_rows_l_kernel(pinned p_dense_matrix<double>& a, const size_t& i_mark, const size_t& k)
{
    scope_select("1 from ambient as remove_rows where master is 0");
    if(!scope.involved()) return;
//    zout << "2dbcd in remove_rows ("<< ambient::rank() <<"):\n"; info(a);

    block_2d_cycle_assign(a);
}

void remove_cols_l_kernel(pinned p_dense_matrix<double>& a, const size_t& j_mark, const size_t& k)
{
    scope_select("1 from ambient as remove_cols where master is 0");
    if(!scope.involved()) return;
//    zout << "2dbcd in remove_cols ("<< ambient::rank() <<"):\n"; info(a);

    block_2d_cycle_assign(a);
}

void sqrt_diagonal_l_kernel(pinned p_dense_matrix<double>& a)
{
    scope_select("1 from ambient as sqrt_diagonal where master is 0");
    if(!scope.involved()) return;
//    zout << "2dbcd in sqrt_diagonal ("<< ambient::rank() <<"):\n"; info(a);

    block_2d_cycle_assign(a);
}


void one_l_scalapack_kernel(const p_dense_matrix<double>& a)
{
    scope_select("1 from ambient as one_scalapack where master is 0");
    if(!scope.involved()) return;
//    zout << "SCALAPACK: 2dbcd in one_scalapack ("<< ambient::rank() <<"):\n"; info(a);

    block_2d_cycle_assign(a);
}

void two_l_scalapack_kernel(const p_dense_matrix<double>& a, const p_dense_matrix<double>& b)
{
    scope_select("1 from ambient as two_scalapack where master is 0");
    if(!scope.involved()) return;
//    zout << "SCALAPACK: 2dbcd in two_scalapack ("<< ambient::rank() <<"):\n"; info(a); info(b);

    block_2d_cycle_assign(a);
    block_2d_cycle_assign(b);
}

void gemm_l_scalapack_kernel(const p_dense_matrix<double>& a, const p_dense_matrix<double>& b,  p_dense_matrix<double>& c)
{
    int num = 1;//get_grid_dim(a).y; 
    scope_select(num+" from ambient as gemm_scalapack where master is 0 and breakdown contains "+get_id(a));
    if(!scope.involved()) return;
//    zout << "2dbcd in gemm_scalapack ("<< ambient::rank() <<"):\n"; info(a); info(b); info(c);

    block_2d_cycle_assign(a);
    block_2d_cycle_assign(b);
    block_2d_cycle_assign(c);
}

void mem_bound_l_kernel(const p_dense_matrix<double>& a, const p_dense_matrix<double>& b,  pinned p_dense_matrix<double>& c)
{
    scope_select(1 +" from ambient as mem_bound where master is 0 and breakdown contains "+ get_id(c));
    if(!scope.involved()) return;
//    zout << "2dbcd in membound ("<< ambient::rank() <<"):\n"; info(a); info(b); info(c);

    block_2d_cycle_assign(a);
    block_2d_cycle_assign(b);
    block_2d_cycle_assign(c);
}

void scale_l_kernel(const p_dense_matrix<double>& m, const double& t, pinned p_dense_matrix<double>& out)
{
    scope_select(1 +" from ambient as scale where master is 0 and breakdown contains "+ get_id(m));
    if(!scope.involved()) return;
//    zout << "2dbcd in scale ("<< ambient::rank() <<"):\n"; info(m); info(out);

    block_2d_cycle_assign(m);
    block_2d_cycle_assign(out);
}
/////////////////////
// testing kernels // 

void svd_l_scalapack_kernel(const p_dense_matrix<double>& a, p_dense_matrix<double>& u, p_dense_matrix<double>& v, p_dense_matrix<double>& s)
{
    int num = 1;
    scope_select(num+" from ambient as svd where master is 0 and breakdown contains "+ get_id(a));
    if(!scope.involved()) return;
//    zout << "2dbcd in svd ("<< ambient::rank() <<"):\n"; info(a); info(u); info(v); info(s);

    block_outright_assign(s);
    block_2d_cycle_assign(a);
    block_2d_cycle_assign(u);
    block_2d_cycle_assign(v);
}

void syev_l_scalapack_kernel(const p_dense_matrix<double>& m, p_dense_matrix<double>& w, p_dense_matrix<double>& z)
{
    int num = 1;
    scope_select(num+" from ambient as syev where master is 0 and breakdown contains "+ get_id(m)); // todo: correct the naming issue
    if(!scope.involved()) return;
//    zout << "2dbcd in syev ("<< ambient::rank() <<"):\n"; info(m); info(w); info(z);

    block_outright_assign(z);
    block_2d_cycle_assign(m);
    block_2d_cycle_assign(w);
}

void gemm_diagonal_lhs_l_kernel(const p_dense_matrix<double>& a_diag, pinned const p_dense_matrix<double>& b, p_dense_matrix<double>& c)
{
    scope_select("1 from ambient as gemm_lhs_diagonal where master is 0");
    if(!scope.involved()) return;
//    zout << "2dbcd in gemm_diagonal_lhs ("<< ambient::rank() <<"):\n"; info(a_diag); info(b); info(c);

    block_2d_cycle_assign(a_diag);
    block_2d_cycle_assign(b);
    block_2d_cycle_assign(c);
}

void gemm_diagonal_rhs_l_kernel(pinned const p_dense_matrix<double>& a, const p_dense_matrix<double>& b_diag, p_dense_matrix<double>& c)
{
    scope_select("1 from ambient as gemm_rhs_diagonal where master is 0");
    if(!scope.involved()) return;
//    zout << "2dbcd in gemm_diagonal_rhs ("<< ambient::rank() <<"):\n"; info(a); info(b_diag); info(c);

    block_2d_cycle_assign(a);
    block_2d_cycle_assign(b_diag);
    block_2d_cycle_assign(c);
}

void validation_l_kernel(pinned const p_dense_matrix<double>& a_ambient, const p_dense_matrix<double>& b_scalapack) 
{ 
    int num = 1; //get_grid_dim(a_ambient).y; 
    scope_select(num+" from ambient as validation where master is 0");
    if(!scope.involved()) return;
//    zout << "2dbcd in validation ("<< ambient::rank() <<"):\n"; info(a_ambient); info(b_scalapack);

    block_2d_cycle_assign(a_ambient); 
    block_2d_cycle_assign(b_scalapack); 
}

void reshape_l2r_l_kernel(const p_dense_matrix<double>& left, pinned p_dense_matrix<double>& right,
                          const size_t& left_offset, const size_t& right_offset, 
                          const size_t& sdim, const size_t& ldim, const size_t& rdim)
{
    int num = 1; //get_grid_dim(a_ambient).y; 
    scope_select(num+" from ambient as reshape_l2r where master is 0 and breakdown contains "+ get_id(left));
    if(!scope.involved()) return;
//    zout << "2dbcd in reshape_l2r ("<< ambient::rank() <<"):\n"; info(left); info(right);

    block_2d_cycle_assign(left); 
    block_2d_cycle_assign(right); 
}

void reshape_r2l_l_kernel(pinned p_dense_matrix<double>& left, const p_dense_matrix<double>& right,
                          const size_t& left_offset, const size_t& right_offset, 
                          const size_t& sdim, const size_t& ldim, const size_t& rdim)
{
    int num = 1; //get_grid_dim(a_ambient).y; 
    scope_select(num+" from ambient as reshape_l2r where master is 0 and breakdown contains "+ get_id(right));
    if(!scope.involved()) return;
//    zout << "2dbcd in reshape_r2l ("<< ambient::rank() <<"):\n"; info(left); info(right);

    block_2d_cycle_assign(left); 
    block_2d_cycle_assign(right); 
}

void associated_validation_l_kernel(pinned const p_dense_matrix<double>& a, const p_dense_matrix<double>& b)
{
    int num = 1;
    scope_select(num+" from ambient as associated_validation where master is 0 and breakdown contains "+ get_id(a)); 
    if(!scope.involved()) return;
    zout << "2dbcd in associated_validation ("<< ambient::rank() <<"):\n"; info(a); info(b);

    block_outright_assign(a);
    block_outright_assign(b);
}

void associated_copy_l_kernel(p_dense_matrix<double>& ac, pinned const p_dense_matrix<double>& a)
{
    int num = 1;
    scope_select(num+" from ambient as associated_copy where master is 0 and breakdown contains "+ get_id(a)); 
    if(!scope.involved()) return;
    zout << "2dbcd in associated_copy ("<< ambient::rank() <<"):\n"; info(ac); info(a);

    block_outright_assign(ac);
    block_outright_assign(a);
}

void associated_sort_l_kernel(pinned p_dense_matrix<double>& a) 
{
    int num = 1;
    scope_select(num+" from ambient as associated_sort where master is 0 and breakdown contains "+ get_id(a));
    if(!scope.involved()) return;
    zout << "2dbcd in assocaited_copy ("<< ambient::rank() <<"):\n"; info(a);

    block_outright_assign(a);
}

void associated_reverse_l_kernel(p_dense_matrix<double>& a, const size_t& num_rows)
{
    int num = 1;
    scope_select(num+" from ambient as associated_reverse where master is 0 and breakdown contains "+ get_id(a));
    if(!scope.involved()) return;
    zout << "2dbcd in associated_reverse ("<< ambient::rank() <<"):\n"; info(a);

    block_outright_assign(a);
}

void touch_l_kernel(p_dense_matrix<double>& target){
    int num = 1; //get_grid_dim(a_ambient).y; 
    scope_select(num+" from ambient as touch_l where master is 0");
    if(!scope.involved()) return;
//    zout << "2dbcd in touch ("<< ambient::rank() <<"):\n"; info(target);

    block_2d_cycle_assign(target); 
}
