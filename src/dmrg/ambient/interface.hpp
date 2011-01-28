#ifndef AMBIENT_INTERFACE_H
#define AMBIENT_INTERFACE_H
#include <assert.h>
#include "ambient/interface/p_profile.h"

namespace ambient {
    scheduler& instance(){ return scheduler::instance(); }

    class charge { // workload of individual rank in terms of workgroups 
    public:
        charge():accept(false){ }
        charge& operator()(const int rank)
        {
            if(rank == 0){ accept = true; }
            else{ accept = false; target = rank; }
            return *this;
        }
        charge& operator+=(workgroup* group)
        {
            if(accept){
                accept = false;
                zout << "I've accepted group" << group->i << group->j << std::endl;
                recvlist.push_back(group);
                group->owner = 0;
            }else if(group->owner == 0){
                group->owner = target;
                sendlist.push_back(group);
            }
        }
        bool accept;
        int target;
        std::list<workgroup*> sendlist;
        std::list<workgroup*> recvlist;
    } charge;
}


namespace blas {

    using namespace ambient;

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

}

namespace ambient {


    using namespace blas;

    template <typename T> p_profile::p_profile(const T* ptr){ 
        p_profile_model(this, ptr); 
        for(int j = 0; j < this->dim.x / scheduler::instance().group_dim().x; j++) // fortran matrix style ,)
            for(int i = 0; i < this->dim.y / scheduler::instance().group_dim().y; i++)
                for(int k = 0; k < this->dim.z / scheduler::instance().group_dim().z; k++)
                    skeleton.push_back(new workgroup(this, i, j, k));
    };

    template <typename RT, typename FC, typename FL, class T1, class T2> 
    const RT push(FC c_kernel, FL l_kernel, const T1& arg1, const T2& arg2){
        zout << "Calling the template that takes 2 args\n";
        p_profile* handle = new p_profile((const RT*)NULL);
        RT out(handle);
        push(c_kernel, l_kernel, arg1, arg2, out);
        return RT(handle);
    }

    template <typename FC, typename FL, class T1, class T2, class T3>
    void push(FC c_kernel, FL l_kernel, const T1& arg1, const T2& arg2, const T3& arg3){
        l_kernel(get_profile(arg1), get_profile(arg2), get_profile(arg3)); 
        c_kernel(get_profile(arg1), get_profile(arg2), get_profile(arg3));
    }
    

    template <typename L, typename R>
    void pin(const L& lhs, const R& rhs){
//        new p_action('=', &lhs, &rhs);
    }

    template <typename L, typename R>
    void copy(const L& lhs, const R& rhs){
//        new p_action('|', &lhs, &rhs);
    }
}

#endif
