#include "ambient/interface/p_profile.h"

namespace ambient {

    workgroup::workgroup(p_profile* p, int i, int j): profile(p), i(i), j(j) {};
    void* workgroup::item(int i, int j = 0, int k = 0){ 
        return NULL;  // this->data[some calculations];
    }

    p_profile::group(int i, int j, int k){
        return skeleton[0*10+0];
    }
