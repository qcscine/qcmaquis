// nested inside ambient.hpp in ambient namespace

template <typename T>  
void breakdown_model(void_pt* profile, const p_dense_matrix<T>* ptr)
{        
    if(ptr == NULL){
        profile->state = PROXY;
        profile->dim.x = 0;
        profile->dim.y = 0;
        profile->dim.z = 0;
    }else{
        profile->init_fp   = matrix_i_kernel;
        profile->t_size    = sizeof(T);
        profile->dim.x     = ptr->num_columns();
        profile->dim.y     = ptr->num_rows();
        profile->dim.z     = 1;
    }
    profile->regroup();
}

template <typename T>
void breakdown_proxy_model(void_pt* proxy, void_pt* profile, const p_dense_matrix<T>* ptr)
{
        proxy->profile      = profile;
        proxy->group_id     = profile->group_id;
        proxy->id           = profile->id;
        proxy->scope        = profile->scope;
        proxy->layout       = profile->layout;         // pointer
        proxy->dim          = profile->dim;
        proxy->t_size       = profile->t_size; 
        proxy->packet_type  = profile->packet_type;    // pointer
        proxy->state        = GENERIC;                 // spike for now
        proxy->imitate(profile); 
        proxy->state        = PROXY;
}

template<>
void_pt& breakdown(const int& obj){ return *(new void_pt(&obj)); }
void breakdown_model(void_pt* profile, const int* ptr){  }//assert(false); }

template<>
void_pt& breakdown(const double& obj){ return *(new void_pt(&obj)); }
void breakdown_model(void_pt* profile, const double* ptr){ assert(false); }

