#include "ambient/ambient.h"
#include "ambient/core/scope_proxy.h"

namespace ambient {

    scope_proxy& scope_proxy::instance()
    {
        static scope_proxy* singleton = NULL;
        if(!singleton) singleton = new scope_proxy();
        return *singleton;
    }
    scope_proxy::scope_proxy(){ }

    scope_proxy& scope_proxy::operator()(const int rank)
    {
        return *this;
    }
    void scope_proxy::set_group(groups::group* grp)
    {
        this->grp = grp;
    }
    groups::group* scope_proxy::get_group(){
        if(this->grp == NULL){
            printf("Attempting to access NULL scope, check if select() was called\n");
        }
        return this->grp;
    }
    int scope_proxy::get_master_g(){
        return this->grp->get_master_g();
    }
    groups::packet_manager* scope_proxy::get_manager(){
        return this->grp->get_manager();
    }
    int scope_proxy::get_size(){
        return this->grp->get_size();
    }
    int scope_proxy::get_rank(){
        return this->grp->get_rank();
    }
    const char* scope_proxy::get_name(){
        return this->grp->get_name();
    }
    void scope_proxy::set_op(core::operation* op){
        this->op = op;
    }
    core::operation* scope_proxy::get_op(){
        return this->op;
    }
    bool scope_proxy::involved(){
        return this->grp->involved();
    }
    bool scope_proxy::is_master(){
        return this->grp->is_master();
    }
}
