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
void_pt& breakdown(const T& obj){
    void_pt** profile_ptr = const_cast<void_pt**>(&obj.profile);
    return *(*profile_ptr = (void_pt*)(*profile_ptr)->dereference());
}

