#include "ambient/ambient.h"
#include "ambient/core/smp.h"

namespace ambient {

    smp& smp::instance()
    {
        static smp* singleton = NULL;
        if(!singleton) singleton = new smp();
        return *singleton;
    }
    smp::smp():interrupt(false){ }

    smp& smp::operator()(const int rank)
    {
        return *this;
    }
    void smp::set_scope(groups::group* scope)
    {
        this->scope = scope;
        if(scope == NULL) return;
        this->scope_size = scope->count;
        this->rank = ambient::rank(this->scope);
    }
    void smp::set_scope(const char* scope)
    {
        if(scope == NULL) this->set_scope((groups::group*)NULL);
        else this->set_scope(groups::group::group_map(scope));
    }
    groups::group* smp::get_scope(){
        if(this->scope == NULL){
            printf("Attempting to access NULL scope, check if select() was called\n");
            throw core::out_of_scope_e();
        }
        return this->scope;
    }
    void smp::trigger_interrupt()
    {
        this->interrupt = !this->interrupt;
    }
}
