#include "ambient/core/operation/operation.pp.hpp"

namespace USER_NAMESPACE{ using namespace ambient;
    template <typename T> 
    void_pt& breakdown(const T& obj){ return breakdown(&obj); }

    template <typename T> 
    void_pt& breakdown(const T* obj){ 
        void_pt** profile_ptr = const_cast<void_pt**>(&obj->profile);
        return *(*profile_ptr = (void_pt*)(*profile_ptr)->dereference());
    }
}

namespace ambient { using namespace USER_NAMESPACE;

    template <typename T> void_pt::void_pt(const T* ptr) : p_profile(){
        breakdown_model(this, ptr);
        this->regroup(); 
    };
}
