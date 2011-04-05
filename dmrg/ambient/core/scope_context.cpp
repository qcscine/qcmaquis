#include "ambient/ambient.h"
#include "ambient/core/scope_context.h"

namespace ambient {

    scope_context& scope_context::instance()
    {
        static scope_context* singleton = NULL;
        if(!singleton) singleton = new scope_context();
        return *singleton;
    }
    scope_context::scope_context():grp(NULL){ }

    scope_context& scope_context::operator()(const int rank)
    {
        return *this;
    }
    void scope_context::set_group(groups::group* grp)
    {
        this->grp = grp;
    }
    void scope_context::reset_group()
    {
        this->grp = groups::group_map("ambient");
    }
    groups::group* scope_context::get_group(){
        if(this->grp == NULL){
            //printf("Attempting to access NULL scope, check if select() was called\n");
        }
        return this->grp;
    }
    int scope_context::get_master_g(){
        return this->grp->get_master_g();
    }
    groups::packet_manager* scope_context::get_manager(){
        return this->grp->get_manager();
    }
    int scope_context::get_size(){
        return this->grp->get_size();
    }
    int scope_context::get_rank(){
        return this->grp->get_rank();
    }
    const char* scope_context::get_name(){
        return this->grp->get_name();
    }
    void scope_context::set_op(core::operation* op){
        this->op = op;
    }
    core::operation* scope_context::get_op(){
        return this->op;
    }
    bool scope_context::involved(){
        return this->grp->involved();
    }
    bool scope_context::is_master(){
        return this->grp->is_master();
    }
}
