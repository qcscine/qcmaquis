#ifndef AMBIENT_INTERFACE_H
#define AMBIENT_INTERFACE_H
#include <assert.h>
#include "ambient/interface/p_profile.h"


namespace blas {
    template <typename T> 
    class p_dense_matrix; // forward declaration of p_dense_matrix

    template <typename T> 
    ambient::p_profile* get_profile(const T* obj){ return obj->profile;      }
    template <typename T>  
    void p_profile_model(ambient::p_profile* profile, const p_dense_matrix<T>* ptr)
    {
        profile->type = "matrix";
        if(ptr == NULL) 
            profile->proxy = true;
        else{
            profile->proxy = false;

            profile->scope = new T[(size_t)(ptr->get_lda()*ptr->get_sda()) + ambient::get_bound()];
            profile->data = (void*)((size_t)profile->scope + ambient::get_bound());

            // TODO: additional construction configuration goes below
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
    }

    void assign_c_kernel(ambient::p_profile* a, ambient::p_profile* b){
        zout << "Executing assign computation kernel...\n";
    }
    void assign_l_kernel(ambient::p_profile* a, ambient::p_profile* b){
        zout << "Executing assign logistics kernel...\n";
    }

}

namespace ambient {

    scheduler& instance(){ return scheduler::instance(); }

    using namespace blas;


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



    template <typename T> p_profile::p_profile(const T* ptr){ 
//        skeleton.push_back(new workgroup(this, 0, 0));
        p_profile_model(this, ptr); 
    };

    template <typename L, typename R>
    p_action::p_action(char op_code, const L* lhs, const R* rhs): arguments(get_profile(lhs), get_profile(rhs)), op_code(op_code) 
    { scheduler::instance().push(this); };


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
        new p_action('+', &arg1, &arg2);
        l_kernel(get_profile(arg1), get_profile(arg2), get_profile(arg3)); 
        c_kernel(get_profile(arg1), get_profile(arg2), get_profile(arg3));
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
