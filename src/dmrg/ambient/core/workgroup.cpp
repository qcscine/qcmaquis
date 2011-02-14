#include "ambient/ambient.h"
#include "ambient/core/workgroup.h"

namespace ambient {

    workgroup::workgroup(p_profile** p, int i, int j, int k): profile(p), i(i), j(j), k(k), initialized(false) {};
    void* workgroup::item(int i, int j, int k){
        p_profile* profile = (*this->profile)->dereference();
        i = i*profile->item_dim().y;
        j = j*profile->item_dim().x;
        k = k*profile->item_dim().z;
    
        int x_size = profile->group_dim().x;
        int y_size = profile->group_dim().y;
        int z_size = profile->group_dim().z;
        if(i >= y_size || j >= x_size || k >= z_size) printf("Warning: accessing group item that is out of range\n");
        return (void*)((size_t)this->data + j*y_size*z_size + i*z_size + k); // the model of accessing actual data can be changed in future - will need to try out!
    }
}
