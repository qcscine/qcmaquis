#include "ambient/ambient.h"
#include "ambient/controllers/workgroup_context.h"

namespace ambient { namespace controllers {

    workgroup_context::workgroup_context(){ 
    }

    void workgroup_context::bind(models::operation* op){
        this->bind_op = op;
       // this->pin = op->pin;
    }

    void workgroup_context::discharge(models::operation* kernel){
        //std::vector<models::imodel::layout::entry>& workload = this->pin->layout->get_list();
        //for(int k=0; k < workload.size(); k++){
        //    this->pin->set_default_block(workload[k].i, workload[k].j);
            kernel->invoke();
            //ambient::spin();
        //}
    }

} }
