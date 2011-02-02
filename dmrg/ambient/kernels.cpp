void plus_c_kernel(ambient::p_profile* a, ambient::p_profile* b, ambient::p_profile* out){
    zout << "Executing plus computation kernel...\n";
}

void info(ambient::p_profile* p){
    zout << "Matrix "<< p->dereference()->id <<" ("<< p->id <<"): size of the task is "<< p->dereference()->grid_dim().y << " x "<< p->dereference()->grid_dim().x << 
            " groups sized "<<p->dereference()->group_dim().y <<" x "<< p->dereference()->group_dim().x  <<
            " items of "<< p->dereference()->item_dim().x << " x " << p->dereference()->item_dim().y << " elements\n";
}

void plus_l_kernel(ambient::p_profile* a, ambient::p_profile* b, ambient::p_profile* out){
    a >> dim3(10,5), dim3(1,1), dim3(10,1);
    info(a); info(b); info(out);

//    select(ALL, ambient);
    charge(0) += a->group(0, 0);
    charge(1) += a->group(0, 1);
    charge(1) += a->group(0, 2);
    charge(1) += a->group(2, 3);
    charge(1) += a->group(3, 4);
}

void assign_c_kernel(ambient::p_profile* a, ambient::p_profile* b){
    zout << "Executing assign computation kernel...\n";
}
void assign_l_kernel(ambient::p_profile* a, ambient::p_profile* b){
    zout << "Executing assign logistics kernel...\n";
}

