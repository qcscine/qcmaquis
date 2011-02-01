void plus_c_kernel(ambient::p_profile* a, ambient::p_profile* b, ambient::p_profile* out){
    zout << "Executing plus computation kernel...\n";
}
void plus_l_kernel(ambient::p_profile* a, ambient::p_profile* b, ambient::p_profile* out){
    a >> dim3(10,5), dim3(1,1), dim3(10,1);
    zout << "Executing plus logistics kernel... size of the task is "<< a->grid_dim().y << " x "<< a->grid_dim().x <<" groups sized "<<a->group_dim().y <<" x "<< a->group_dim().x  <<" items of "<< a->item_dim().x << " x " << a->item_dim().y << " elements\n";
    zout << "Executing plus logistics kernel... size of the task is "<< b->grid_dim().y << " x "<< b->grid_dim().x <<" groups sized "<<b->group_dim().y <<" x "<< b->group_dim().x  <<" items of "<< b->item_dim().x << " x " << b->item_dim().y << " elements\n";
    zout << "Executing plus logistics kernel... size of the task is "<< out->grid_dim().y << " x "<< out->grid_dim().x <<" groups sized "<<out->group_dim().y <<" x "<< out->group_dim().x  <<" items of "<< out->item_dim().x << " x " << out->item_dim().y << " elements\n";
    charge(0) += a->group(0, 0);
    charge(0) += a->group(0, 1);
}

void assign_c_kernel(ambient::p_profile* a, ambient::p_profile* b){
    zout << "Executing assign computation kernel...\n";
}
void assign_l_kernel(ambient::p_profile* a, ambient::p_profile* b){
    zout << "Executing assign logistics kernel...\n";
}

