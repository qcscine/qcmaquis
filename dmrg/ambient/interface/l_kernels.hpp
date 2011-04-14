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
        printf("Matrix %d:%d size of the task is %d x %d groups sized %d x %d items of %d x %d elements\n", 
               *p.group_id, p.id, p.get_grid_dim().y, p.get_grid_dim().x, p.get_group_dim().y, p.get_group_dim().x, p.get_item_dim().x, p.get_item_dim().y);
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

void check_gemm_l_kernel(pinned const p_dense_matrix<double>& c)
{
    scope_select("2 from ambient as gemm_check where master is 0 and breakdown contains "+get_id(c)); // todo: correct the naming issue
    if(!scope.involved()) return;

    zout << "2d-block-cyclic decomposition kernel in gemm ("<< ambient::rank() <<"):\n"; info(c);
    block_2d_cycle_assign(c);
}

void gemm_l_kernel(pinned const p_dense_matrix<double>& a, const p_dense_matrix<double>& b, p_dense_matrix<double>& c)
{
    //breakdown(a) >> dim3(), dim3(), dim3();
    scope_select("2 from ambient as gemm where master is 0"); // todo: correct the naming issue
    if(!scope.involved()) return;

    zout << "2d-block-cyclic decomposition kernel in gemm ("<< ambient::rank() <<"):\n"; info(a); info(b); info(c);

    block_2d_cycle_assign(a);
    block_2d_cycle_assign(b);
    block_2d_cycle_assign(c);
}

void copy_l_kernel(p_dense_matrix<double>& ac, pinned const p_dense_matrix<double>& a)
{
    scope_select("2 from ambient as copy where master is 0");
    if(!scope.involved()) return; // out of scope quick exit
    zout << "2d-block-cyclic decomposition kernel in copy ("<< ambient::rank() <<"):\n"; info(a); info(ac);

    block_2d_cycle_assign(a);
    block_2d_cycle_assign(ac);
}

void copy_l_kernel3(p_dense_matrix<double>& ac, pinned const p_dense_matrix<double>& a, p_dense_matrix<double>& d)
{
    scope_select("2 from ambient as copy where master is 0");
    if(!scope.involved()) return; // out of scope quick exit
    zout << "2d-block-cyclic decomposition kernel in copy ("<< ambient::rank() <<"):\n"; info(a); info(ac); info(d);

    block_2d_cycle_assign(a);
    block_2d_cycle_assign(ac);
    block_2d_cycle_assign(d);
}

void resize_l_kernel(p_dense_matrix<double>& a, const size_t& rows, const size_t& cols)
{
    breakdown(a).set_dim(ambient::dim3(cols,rows));

    scope_select("2 from ambient as resize where master is 0");
    if(!scope.involved()) return; // out of scope quick exit

    zout << "2d-block-cyclic decomposition kernel in resize matrix ("<< ambient::rank() <<"):\n"; info(a);
    block_2d_cycle_assign(a);
}

void remove_rows_l_kernel(pinned p_dense_matrix<double>& a, const size_t& i_mark, const size_t& k)
{
    scope_select("2 from ambient as remove_rows where master is 0");
    if(!scope.involved()) return; // out of scope quick exit
    zout << "2d-block-cyclic decomposition kernel in remove rows ("<< ambient::rank() <<"):\n"; info(a);
    block_2d_cycle_assign(a);
}

void remove_cols_l_kernel(pinned p_dense_matrix<double>& a, const size_t& j_mark, const size_t& k)
{
    scope_select("2 from ambient as remove_cols where master is 0");
    if(!scope.involved()) return; // out of scope quick exit
    zout << "2d-block-cyclic decomposition kernel in remove cols ("<< ambient::rank() <<"):\n"; info(a);
    block_2d_cycle_assign(a);
}

void sqrt_diagonal_l_kernel(pinned p_dense_matrix<double>& a)
{
    scope_select("2 from ambient as sqrt_diagonal where master is 0");
    if(!scope.involved()) return; // out of scope quick exit
    zout << "2d-block-cyclic decomposition kernel in sqrt diagonal  ("<< ambient::rank() <<"):\n"; info(a);
    block_2d_cycle_assign(a);
}


void one_l_scalapack_kernel(const p_dense_matrix<double>& a)
{
    scope_select("* from ambient as one_scalapack where master is 0");
    scope_retain("* from ambient as work_storage");
    if(!scope.involved()) return; // out of scope quick exit
    zout << "2d-block-cyclic decomposition kernel in one scalapack kernel ("<< ambient::rank() <<"):\n"; info(a);
    block_2d_cycle_assign(a);
}

void two_l_scalapack_kernel(const p_dense_matrix<double>& a, const p_dense_matrix<double>& b)
{
    scope_select("* from ambient as two_scalapack where master is 0");
    scope_retain("* from ambient as work_storage");
    if(!scope.involved()) return; // out of scope quick exit
    block_2d_cycle_assign(a);
    block_2d_cycle_assign(b);
}

void three_l_scalapack_kernel(const p_dense_matrix<double>& a, const p_dense_matrix<double>& b,  const p_dense_matrix<double>& c)
{
    scope_select("* from ambient as three_scalapack where master is 0");
    scope_retain("* from ambient as work_storage");
    if(!scope.involved()) return; // out of scope quick exit
    block_2d_cycle_assign(a);
    block_2d_cycle_assign(b);
    block_2d_cycle_assign(c);
}


void gemm_l_scalapack_kernel(const p_dense_matrix<double>& a, const p_dense_matrix<double>& b,  p_dense_matrix<double>& c)
{
    scope_select("* from ambient as gemm_scalapack where master is 0");
    scope_retain("* from ambient as work_storage");
    if(!scope.involved()) return; // out of scope quick exit
    zout << "2d-block-cyclic decomposition kernel in gemm_scalapack ("<< ambient::rank() <<"):\n"; info(a); info(b); info(c);

    block_2d_cycle_assign(a);
    block_2d_cycle_assign(b);
    block_2d_cycle_assign(c);
}

void mem_bound_l_kernel(const p_dense_matrix<double>& a, const p_dense_matrix<double>& b,  pinned p_dense_matrix<double>& c)
{
    scope_select(3 +" from ambient as mem_bound where master is "+ 1 +" and breakdown contains "+ get_id(c));
    scope_retain("2 from ambient as work_storage");
    if(!scope.involved()) return; // out of scope quick exit
    zout << "2d-block-cyclic decomposition kernel in membound ("<< ambient::rank() <<"):\n"; info(a); info(b); info(c);

    block_2d_cycle_assign(a);
    block_2d_cycle_assign(b);
    block_2d_cycle_assign(c);
}

void scale_l_kernel(const p_dense_matrix<double>& m, const double& t, pinned p_dense_matrix<double>& out)
{
    scope_select(3 +" from ambient as scale where master is "+ 1 +" and breakdown contains "+ get_id(m));
    scope_retain("2 from ambient as work_storage");
    if(!scope.involved()) return; // out of scope quick exit
    zout << "2d-block-cyclic decomposition kernel in scale ("<< ambient::rank() <<"):\n"; info(m); info(out);

    block_2d_cycle_assign(m);
    block_2d_cycle_assign(out);
}
/////////////////////
// testing kernels // 

//universal logistic kernel, distribution should be donne before

void null_l_scalapack_svd_kernel(const p_dense_matrix<double>  &  M, p_dense_matrix<double>  & U, p_dense_matrix<double> & V,  double* & S )
{
    scope_select("* from ambient as null_l_scalapack_svd where master is 0");
    scope_retain("* from ambient as work_storage");
    if(!scope.involved()) return; // out of scope quick exit

    block_2d_cycle_assign(M);
    block_2d_cycle_assign(U);
    block_2d_cycle_assign(V);

}

/*
void null_l_scalapack_svd_kernel(const p_dense_matrix<double>  &  M, p_dense_matrix<double>  & U, p_dense_matrix<double> & V,  p_dense_matrix<double> & S )
{
    scope_select("* from ambient as work where master is 0");
    scope_retain("* from ambient as work_storage");
   
    if(!scope.involved()) return; // out of scope quick exit

    block_2d_cycle_assign(M);
    block_2d_cycle_assign(U);
    block_2d_cycle_assign(V);

}// do nothing just need to test
*/




// pdiag*pdense
void gemm_lhs_diagonal_l_kernel(const p_dense_matrix<double> & a_diag, pinned const p_dense_matrix<double>& b,  p_dense_matrix<double>&  c)
{
    scope_select("* from ambient as gemm_lhs_diagonal where master is 0");
    scope_retain("* from ambient as work_storage");
    if(!scope.involved()) return; // out of scope quick exit

    block_2d_cycle_assign(a_diag);
    block_2d_cycle_assign(b);
    block_2d_cycle_assign(c);
}

// pdense*pdiag
void gemm_rhs_diagonal_l_kernel(pinned const p_dense_matrix<double> & a, const p_dense_matrix<double>& b_diag,  p_dense_matrix<double>&  c)
{
    scope_select("* from ambient as gemm_rhs_diagonal where master is 0");
    scope_retain("* from ambient as work_storage");
    if(!scope.involved()) return; // out of scope quick exit

    block_2d_cycle_assign(a);
    block_2d_cycle_assign(b_diag);
    block_2d_cycle_assign(c);
}

void copy_svd_l_kernel(pinned p_dense_matrix<double> & a, double* & Array) 
{ 
    scope_select("* from ambient as copy_svd where master is 0"); 
    scope_retain("* from ambient as work_storage"); 
    if(!scope.involved()) return; // out of scope quick exit 
 
    block_2d_cycle_assign(a); 
} 

void init_double_l_kernel(pinned p_dense_matrix<double> & a)
{
    scope_select("* from ambient as init_double where master is 0");
    scope_retain("* from ambient as work_storage");
    if(!scope.involved()) return; // out of scope quick exit

    block_2d_cycle_assign(a);
}

