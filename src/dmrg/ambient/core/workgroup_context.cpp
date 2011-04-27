#include "ambient/ambient.h"
#include "ambient/core/workgroup_context.h"

namespace ambient {

    workgroup_context::workgroup_context(){ }

    void workgroup_context::bind(p_profile* pin, int i, int j){
        this->pin = pin;
        this->i = i;
        this->j = j;
    }

    void discharge(core::operation* kernel){

    }

}
