#include "ambient/ambient.h"
#include "ambient/core/smp.h"

namespace ambient {

typedef ambient::p_profile   void_pt;
typedef ambient::p_profile_s void_spt;

    smp& smp::instance()
    {
        static smp* singleton = NULL;
        if(!singleton) singleton = new smp();
        return *singleton;
    }
    smp::smp():accept(true){ }

    smp& smp::operator()(const int rank)
    {
        return *this;
    }
    void smp::assign(workgroup* group)
    {
        if(this->accept){
            if(this->id == UNDEFINED_ID) return;
//            printf("%s: p%d: I've accepted group %d %d\n", this->scope->name, this->id, group->i,group->j);
            recvlist.push_back(group);
            group->owner = ambient::rank();
        }else{
            if(group->owner == ambient::rank())
            { // this place doesn't allow to assign the same group to different procs
                group->owner = this->scope->translate_rank(this->id);
                sendlist.push_back(group);
            }
        }
    }
    void smp::set_scope(const groups::group* scope)
    {
        this->scope = scope;
        this->scope_size = scope->count;
        this->id = ambient::rank(this->scope);
    }
    void smp::set_scope(const char* scope)
    {
        this->scope = groups::group::group_map(scope);
        this->scope_size = this->scope->count;
        this->id = ambient::rank(this->scope);
    }
    void smp::get_info(core::operation* op)
    { // to be removed if to support user flexibility in logistics kernels
        int latch = this->id;
        this->accept = false;
        for(int i=0; i < this->scope_size; i++){
            if(i != latch){
                this->id = i;
                op->perform();
            }
        }
        this->id = latch;
        this->accept = true;
    }

    void assign(void_spt* ptr, int i, int j, int k)
    {
        asmp.assign(ptr->group(i, j, k));
    }

}
