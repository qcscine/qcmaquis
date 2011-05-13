#include "ambient/ambient.h"
#include "ambient/core/workgroup_context.h"

namespace ambient {

    workgroup_context::workgroup_context(){ }

    void workgroup_context::bind(core::operation* op){
        this->bind_op = op;
        this->pin = op->pin;
    }

    void workgroup_context::discharge(core::operation* kernel){
        std::vector<core::layout_table::entry>& workload = this->pin->layout->get_list();
        for(int k=0; k < workload.size(); k++){
            this->pin->set_default_block(workload[k].i, workload[k].j);
            kernel->invoke();
            ambient::spin();
        }
    }
    void workgroup_context::finalize(){
        this->bind_op->finalize();
    }
}
