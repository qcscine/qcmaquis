// nested inside ambient.hpp in ambient namespace
#include "ambient/core/operation/operation.pp.hpp"
class void_pt: public p_profile 
{ 
public: 
    template <typename T> void_pt(const T* ptr) : p_profile()
    { breakdown_model(this, ptr); }
private: 
    ~void_pt(){ };  // making user unable to delete the profile
    friend class T; // the container can delete its profile
    template<class T> friend inline void boost::checked_delete(T * x); // so as boost >_<
};

template <typename T> 
void_pt& breakdown(T& obj){
    void_pt** profile_ptr = const_cast<void_pt**>(&obj.profile);
    *profile_ptr = (void_pt*)(*profile_ptr)->dereference();
    (*profile_ptr)->inconstant();
    return **profile_ptr;
}

template <typename T> 
void_pt& breakdown(const T& obj){
    void_pt** profile_ptr = const_cast<void_pt**>(&obj.profile);
    *profile_ptr = (void_pt*)(*profile_ptr)->dereference();
    (*profile_ptr)->constant();
    return **profile_ptr;
}

// aliases for computational kernels
template <typename T>
void_pt& current(T& obj){
    return breakdown(obj);
}
template <typename T>
void_pt& reduced(T& obj, char R){
    if(breakdown(obj).associated_proxy == NULL){
        breakdown(obj).associate_proxy(new void_pt((T*)NULL), R);
        breakdown_proxy_model((void_pt*)breakdown(obj).associated_proxy, &(breakdown(obj)), &obj);
    }
    //return breakdown(obj);
    return *((void_pt*)breakdown(obj).associated_proxy);
}
