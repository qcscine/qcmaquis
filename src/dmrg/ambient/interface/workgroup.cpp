#include "ambient/interface/workgroup.h"

namespace ambient {

    workgroup::workgroup(p_profile* p, int i, int j, int k): profile(p), i(i), j(j), k(k) {};
    void* workgroup::item(int i, int j, int k){ 
        i = i*scheduler::instance().item_dim().y;
        j = j*scheduler::instance().item_dim().x;
        k = k*scheduler::instance().item_dim().z;
    
        int x_size = scheduler::instance().group_dim().x;
        int y_size = scheduler::instance().group_dim().y;
        int z_size = scheduler::instance().group_dim().z;
        if(i >= y_size || j >= x_size || k >= z_size) printf("Warning: accessing group item that is out of range\n");
        return (void*)((size_t)this->data + j*y_size*z_size + i*z_size + k); // the model of accessing actual data can be changed in future - will need to try out!
    }
}
