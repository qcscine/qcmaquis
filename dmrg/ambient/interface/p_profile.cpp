#include "ambient/interface/p_profile.h"

namespace ambient {

    workgroup::workgroup(p_profile* p, int i, int j, int k): profile(p), i(i), j(j), k(k) {};
    void* workgroup::item(int i, int j, int k){ 
        return NULL; // this->data[some calculations]; // the model of accessing actual data can be changed in future - will need to try out!
    }

    workgroup* p_profile::group(int i, int j, int k){
        if(this->proxy){
            return new workgroup(this, i, j, k);
        }else{
            int x_size = this->dim.x / scheduler::instance().group_dim().x;
            int y_size = this->dim.y / scheduler::instance().group_dim().y;
            int z_size = this->dim.z / scheduler::instance().group_dim().z;

            if(i > y_size || j > x_size || k > z_size) printf("Warning: accessing group that is our of range\n");
            return skeleton[ j*y_size*z_size + i*z_size + k  ];
        }
    }
}
