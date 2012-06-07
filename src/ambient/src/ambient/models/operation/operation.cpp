#include "ambient/ambient.h"
#include "ambient/models/operation/operation.h"

namespace ambient { namespace models {

    operation::~operation(){
        (this->*cleanup)();
        free(this->arguments);
        this->grp->idle();
    }

    void operation::invoke(){
        (this->*prototype)();
    }

    void operation::weight(){
        (this->*creditup)();
    }

    void operation::set_group(channels::group* grp){
        this->grp = grp;
        (this->*place)();
    }

    channels::group* operation::get_group(){
        return this->grp;
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

    imodel::revision* operation::get_pin(){
        return this->pin;
    }

    bool operation::pretend(){
        bool pretend = false;
        if(this->pin == NULL){
            pthread_mutex_lock(&this->mutex);
            if(--this->workload) pretend = true;
            pthread_mutex_unlock(&this->mutex);
        }
        return pretend;
    }

    void operation::add_condition(){
        pthread_mutex_lock(&this->mutex);
        this->workload++;
        pthread_mutex_unlock(&this->mutex);
    }

} }
