typedef ambient::p_profile   void_pt;
typedef ambient::p_profile_s void_spt;
 
void info(void_pt* p){
    if(!asmp.accept) return;
    zout << "Matrix "<< p->dereference()->id <<" ("<< p->id <<"): size of the task is "<< p->dereference()->grid_dim().y << " x "<< p->dereference()->grid_dim().x << 
            " groups sized "<<p->dereference()->group_dim().y <<" x "<< p->dereference()->group_dim().x  <<
            " items of "<< p->dereference()->item_dim().x << " x " << p->dereference()->item_dim().y << " elements\n";
}

void plus_l_kernel(void_pt* a, void_pt* b, void_spt* out){
    a >> dim3(10,5), dim3(1,1), dim3(10,1);
    info(a); info(b); info((void_pt*)out);

    select("0.5 from ambient as work");

    for(int i=0; i < out->grid_dim().y; i++)
        for(int j=0; j < out->grid_dim().x; j++)
            if(j % asmp.scope_size == asmp.id)
                assign(a, i, j);
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

