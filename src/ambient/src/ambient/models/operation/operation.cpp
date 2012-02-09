#include "ambient/ambient.h"
#include "ambient/models/operation/operation.h"

namespace ambient { namespace models {

    operation::~operation(){
        (this->*cleanup)();
        free(this->arguments);
        this->grp->idle();
    }
    void operation::invoke(){
        if(this->state == MARKUP){ 
            this->state = LOGISTICS;
        }else if(this->state == LOGISTICS){ 
            this->state = COMPUTING; 
            this->op = this->computing_ptr; 
        }

        (this->*prototype)();
    }
    void operation::weight(){
        (this->*creditup)();
    }
    void operation::set_group(channels::group* grp){
        this->grp = grp;
    }
    size_t operation::get_weight(){
        return this->credit;
    }
    void operation::set_weight(size_t credit){
        this->credit = credit;
    }
    imodel::object& operation::get_vellum(){
        return *this->vellum;
    }
    void operation::set_vellum(imodel::object& v){
        this->vellum = &v;
    }
} }
