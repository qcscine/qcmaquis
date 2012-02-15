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
            this->workload = this->pin->get_layout().
                             get_grid_dim().square();
            pthread_mutex_init(&this->mutex, NULL);
        }

        (this->*prototype)();

        if(this->state == COMPUTING){
            pthread_mutex_lock(&this->mutex);
            if(--this->workload == 0) 
                controller.atomic_complete();
            pthread_mutex_unlock(&this->mutex);
        }
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
    imodel::revision& operation::get_vellum(){
        return *this->vellum;
    }
    void operation::set_vellum(imodel::revision& v){
        this->vellum = &v;
    }
    imodel::revision& operation::get_pin(){
        return *this->pin;
    }

} }
