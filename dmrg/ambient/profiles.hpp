template <typename T> 
class p_dense_matrix; // forward declaration of p_dense_matrix

template <typename T> 
ambient::p_profile* get_profile(const T* obj){ return obj->profile;      }
template <typename T>  
void p_profile_model(ambient::p_profile* profile, const p_dense_matrix<T>* ptr)
{
    profile->type = "matrix";
    if(ptr == NULL){
	profile->proxy = true;
	profile->dim.x = 0;
	profile->dim.y = 0;
	profile->dim.z = 0;
    }else{
	profile->proxy = false;
	profile->scope = new T[(size_t)(ptr->get_lda()*ptr->get_sda()) + ambient::get_bound()];
	profile->data = (void*)((size_t)profile->scope + ambient::get_bound());
	profile->dim.x = ptr->num_columns();
	profile->dim.y = ptr->num_rows();
	profile->dim.z = 1;
    }
}

ambient::p_profile* get_profile(const int* obj){ return new ambient::p_profile(obj); }
void p_profile_model(ambient::p_profile* profile, const int* ptr){ assert(false); }

ambient::p_profile* get_profile(const double* obj){ return new ambient::p_profile(obj); }
void p_profile_model(ambient::p_profile* profile, const double* ptr){ assert(false); }

template <typename T> 
ambient::p_profile* get_profile(const T& obj){ return get_profile(&obj); }

