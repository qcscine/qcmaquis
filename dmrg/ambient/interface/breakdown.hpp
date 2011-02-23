namespace USER_NAMESPACE{ using namespace ambient;

    template <typename T> 
    class p_dense_matrix; // forward declaration of p_dense_matrix
    
    void matrix_i_kernel(ambient::workgroup* grp){
        // dumb 0-initialization for the start >_< 
        memset(grp->data, 0, grp->get_group_dim().y*grp->get_item_dim().y*grp->get_group_dim().x*grp->get_item_dim().x*grp->get_profile()->type_size);
    }
    
    template <typename T>  
    void breakdown_model(ambient::void_pt* profile, const p_dense_matrix<T>* ptr)
    {        
        profile->type = "matrix";
        if(ptr == NULL){
            profile->proxy = true;
            profile->scope = new T[ambient::get_bound()];
            profile->dim.x = 0;
            profile->dim.y = 0;
            profile->dim.z = 0;
        }else{
            profile->proxy   = false;
            profile->init_fp = matrix_i_kernel;
            profile->type_size = sizeof(T);
            profile->dim.x   = ptr->num_columns();
            profile->dim.y   = ptr->num_rows();
            profile->dim.z   = 1;
        }
    }
    
    ambient::void_pt& breakdown(const int* obj){ return *(new ambient::void_pt(obj)); }
    void breakdown_model(ambient::void_pt* profile, const int* ptr){  }//assert(false); }

    ambient::void_pt& breakdown(const double* obj){ return *(new ambient::void_pt(obj)); }
    void breakdown_model(ambient::void_pt* profile, const double* ptr){ assert(false); }
    
} 
