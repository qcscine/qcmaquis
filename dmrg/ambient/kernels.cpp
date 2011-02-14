typedef ambient::p_profile   void_pt;
typedef ambient::p_profile_s void_spt;

void* matrix_i_kernel(ambient::workgroup* grp){ 

    return NULL; 

}
 
void info(void_pt* p){
    if(rank.is_master(asmp.get_scope()))
    printf("Matrix %d:%d size of the task is %d x %d groups sized %d x %d items of %d x %d elements\n", *p->group_id, p->id, p->grid_dim().y, p->grid_dim().x, p->group_dim().y, p->group_dim().x, p->item_dim().x, p->item_dim().y);
}

void plus_l_kernel(void_pt* a, void_pt* b, void_spt* out){
//    a >> dim3(10,5), dim3(1,1), dim3(10,1); <- kinda non-trivial - need to think

    select("0.5 from ambient as work where master is 0");
    info(a); info(b); info((void_pt*)out);

    for(int i=0; i < out->grid_dim().y; i++)
        for(int j=0; j < out->grid_dim().x; j++)
            if(j % asmp.scope_size == asmp.rank){
                assign(a, i, j);
                assign(b, i, j);
                assign(out, i, j);
            }
}

void plus_c_kernel(void_pt* a, void_pt* b, void_spt* out){
    zout << "Executing plus computation kernel...\n";
}


void assign_c_kernel(void_pt* a, void_pt* b){
    zout << "Executing assign computation kernel...\n";
}
void assign_l_kernel(void_pt* a, void_pt* b){
    zout << "Executing assign logistics kernel...\n";
}

