namespace blas{ using namespace ambient;

template <typename T> 
class p_dense_matrix; // forward declaration of p_dense_matrix

void matrix_i_kernel(ambient::workgroup* grp){
        // dumb 0-initialization for the start >_< 
        memset(grp->data, 0, grp->get_group_dim().y*grp->get_item_dim().y*grp->get_group_dim().x*grp->get_item_dim().x*grp->get_profile()->type_size);
}

void_pt* dereference(void_pt** profile_ptr)
{
    *profile_ptr = (void_pt*)(*profile_ptr)->dereference();
    return *profile_ptr;
}

template <typename T> 
ambient::void_pt* get_profile(const T* obj){ return dereference(const_cast<ambient::void_pt**>(&obj->profile)); }

template <typename T>  
void void_pt_model(ambient::void_pt* profile, const p_dense_matrix<T>* ptr)
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

ambient::void_pt* get_profile(const int* obj){ return new ambient::void_pt(obj); }
void void_pt_model(ambient::void_pt* profile, const int* ptr){  }//assert(false); }

ambient::void_pt* get_profile(const double* obj){ return new ambient::void_pt(obj); }
void void_pt_model(ambient::void_pt* profile, const double* ptr){ assert(false); }

template <typename T> 
ambient::void_pt* get_profile(const T& obj){ return get_profile(&obj); }

} namespace ambient { using namespace blas;
    template <typename T> void_pt::void_pt(const T* ptr) : p_profile()
    {
        void_pt_model(this, ptr);
        this->regroup(); 
    };
}
