#include "ambient/ambient.h"
#include "ambient/core/workgroup.h"

namespace ambient {

    workgroup::workgroup(p_profile** p, int i, int j, int k)
    : profile(p), i(i), j(j), k(k), header(NULL), data(NULL), timestamp(0) {};

    p_profile* workgroup::get_profile(){
        return *this->profile;
    }

    dim3 workgroup::get_group_dim(){
        return this->get_profile()->get_group_dim();
    }

    dim3 workgroup::get_item_dim(){
        return this->get_profile()->get_item_dim();
    }

    void workgroup::set_memory(void* memory){
        this->header = memory;
        this->data = (void*)((size_t)memory + this->get_profile()->get_bound());
        this->timestamp = this->get_profile()->timestamp;
    }

    void* workgroup::item(int i, int j, int k){
        p_profile* profile = *this->profile;
        i = i*profile->get_item_dim().y;
        j = j*profile->get_item_dim().x;
        k = k*profile->get_item_dim().z;
    
        int x_size = profile->get_group_dim().x;
        int y_size = profile->get_group_dim().y;
        int z_size = profile->get_group_dim().z;
        if(i >= y_size || j >= x_size || k >= z_size) printf("Warning: accessing group item that is out of range (%d %d %d)\n", i, j, k);
        return (void*)((size_t)this->data + j*y_size*z_size + i*z_size + k); // the model of accessing actual data can be changed in future - will need to try out!
    }
}
