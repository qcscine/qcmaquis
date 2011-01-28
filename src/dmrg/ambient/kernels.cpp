void plus_c_kernel(ambient::p_profile* a, ambient::p_profile* b, ambient::p_profile* out){
    zout << "Executing plus computation kernel...\n";
}
void plus_l_kernel(ambient::p_profile* a, ambient::p_profile* b, ambient::p_profile* out){
    zout << "Executing plus logistics kernel...\n";
    charge(0) += a->group(0, 0);
    charge(0) += a->group(0, 1);
}

void assign_c_kernel(ambient::p_profile* a, ambient::p_profile* b){
    zout << "Executing assign computation kernel...\n";
}
void assign_l_kernel(ambient::p_profile* a, ambient::p_profile* b){
    zout << "Executing assign logistics kernel...\n";
}

