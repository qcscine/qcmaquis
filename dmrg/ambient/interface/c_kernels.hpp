void add_c_kernel(const p_dense_matrix<double>& a, const p_dense_matrix<double>& b, pinned p_dense_matrix<double>& out){
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

void sub_c_kernel(const p_dense_matrix<double>& a, const p_dense_matrix<double>& b, pinned p_dense_matrix<double>& out){
// todo
}

void gemm_c_kernel(const p_dense_matrix<double>& a, const p_dense_matrix<double>& b, pinned p_dense_matrix<double>& out){
// todo
}

void scale_c_kernel(const p_dense_matrix<double>& a, const double& b, pinned p_dense_matrix<double>& out){
// todo
}

void null_c_kernel(p_dense_matrix<double>& a){
    printf("R%d: Executing NULL kernel\n", scope.get_rank());
}

void single_integer_c_kernel(int& input){
    input += 13;
    zout << "single integer kernel: output is " << input << "\n";
}
