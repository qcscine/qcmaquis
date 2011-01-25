#ifndef AMBIENT_INTERFACE_H
#define AMBIENT_INTERFACE_H
#include "ambient/interface/p_profile.h"

namespace blas {
    template <typename T> 
    class p_dense_matrix; // forward declaration of p_dense_matrix

    template <typename T> 
    ambient::p_profile* get_profile(const T* obj)
    { 
        return obj->profile; 
    }
    template <typename T>  
    void p_profile_model(ambient::p_profile* profile, const p_dense_matrix<T>* ptr)
    {
        profile->type = "matrix";
        if(ptr == NULL) 
            profile->proxy = true;
        else{
            profile->proxy = false;
            // TODO: additional construction configuration goes below
        }
    }

    ambient::p_profile* get_profile(const int* obj){ return new ambient::p_profile(obj); }
    void p_profile_model(ambient::p_profile* profile, const int* ptr){ printf("construct_p_profile for int\n");  }

    ambient::p_profile* get_profile(const double* obj){ return new ambient::p_profile(obj); }
    void p_profile_model(ambient::p_profile* profile, const double* ptr){ printf("construct_p_profile for double\n");  }
}

namespace ambient {

    using namespace blas;

    template <typename T> p_profile::p_profile(const T* ptr){ p_profile_model(this, ptr); };

    template <typename L, typename R>
    p_action::p_action(char op_code, const L* lhs, const R* rhs): arguments(get_profile(lhs), get_profile(rhs)), op_code(op_code) {
        scheduler::instance().push(this);
    };

    template <typename T, typename L, typename R>
    const T push(char op_code, const L& lhs, const R& rhs){
        p_profile* handle = new p_profile((const T*)NULL);
        handle->action = new p_action(op_code, &lhs, &rhs);
        return T(handle); 
    }

    template <typename L, typename R>
    void pin(const L& lhs, const R& rhs){
        new p_action('=', &lhs, &rhs);
    }

    template <typename L, typename R>
    void copy(const L& lhs, const R& rhs){
        new p_action('|', &lhs, &rhs);
    }
}

#endif
