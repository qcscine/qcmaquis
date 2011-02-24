#include "ambient/ambient.h"
#include "ambient/core/smp.h"

namespace ambient {

    smp& smp::instance()
    {
        static smp* singleton = NULL;
        if(!singleton) singleton = new smp();
        return *singleton;
    }
    smp::smp(){ }

    smp& smp::operator()(const int rank)
    {
        return *this;
    }
    void smp::set_group(groups::group* group)
    {
        this->group = group;
        if(group == NULL){
            this->rank = UNDEFINED_RANK; return;
        }
        this->size = group->count;
        this->rank = ambient::rank(this->group);
    }
    void smp::set_group(const char* group)
    {
        if(group == NULL) this->set_group((groups::group*)NULL);
        else this->set_group(groups::group::group_map(group));
    }
    groups::group* smp::get_group(){
        if(this->group == NULL){
            printf("Attempting to access NULL scope, check if select() was called\n");
            throw core::out_of_scope_e();
        }
        return this->group;
    }
    int smp::get_size(){
        return this->size;
    }
    int smp::get_rank(){
        return this->rank;
    }
    void smp::set_op(core::operation* op){
        this->op = op;
    }
    core::operation* smp::get_op(){
        return this->op;
    }
    bool smp::involved(){
        return ((this->rank != MPI_UNDEFINED) ? true : false); 
    }
}
