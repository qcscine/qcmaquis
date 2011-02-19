#define __a_once__ asmp.trigger_interrupt();

void matrix_i_kernel(ambient::workgroup* grp){
        // dumb 0-initialization for the start >_< 
        memset(grp->data, 0, grp->get_group_dim().y*grp->get_item_dim().y*grp->get_group_dim().x*grp->get_item_dim().x*grp->get_profile()->type_size);
}
 
void info(const void_pt& p){
    if(rank.is_master(asmp.get_scope()))
    printf("Matrix %d:%d size of the task is %d x %d groups sized %d x %d items of %d x %d elements\n", *p.group_id, p.id, p.get_grid_dim().y, p.get_grid_dim().x, p.get_group_dim().y, p.get_group_dim().x, p.get_item_dim().x, p.get_item_dim().y);
}

void plus_l_kernel(const void_pt& a, void_pt& b, void_spt& out){
//    a >> dim3(10,5), dim3(1,1), dim3(10,1); <- kinda non-trivial - need to think
    select("0.5 from ambient as work where master is 0");
    retain("2 from ambient as work_storage");

    info(a); info(b); info((void_pt&)out);

    for(int i=0; i < out.get_grid_dim().y; i++)
        for(int j=0; j < out.get_grid_dim().x; j++)
            if(j % asmp.scope_size == asmp.rank){
                assign(a,   i, j);
                assign(b,   i, j);
                assign(out, i, j);
            }
}

void plus_c_kernel(const void_pt& a, void_pt& b, void_spt& out){
//    __a_once__
    double* output = out;
    double* ad = a(out.get_group_id().x, out.get_group_id().y);
    double* bd = b(out.get_group_id().x, out.get_group_id().y);
    int size = out.get_group_dim().x*out.get_item_dim().x*
               out.get_group_dim().y*out.get_item_dim().y;
//    printf("R%d: Executing plus computation kernel (%d ops)... for out grp %d %d\n", asmp.rank, size, out.get_group_id().x, out.get_group_id().y);
//    for(int i=0; i < size; i++){
//        output[i] = ad[i]+bd[i];
//    }
}


void assign_c_kernel(void_pt& a, void_pt& b){
    zout << "Executing assign computation kernel...\n";
}
void assign_l_kernel(void_pt& a, void_pt& b){
    zout << "Executing assign logistics kernel...\n";
}

