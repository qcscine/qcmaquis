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
    void smp::assign(workgroup* group)
    {
        if(this->rank == UNDEFINED_RANK) return;
//        printf("%s: p%d: I've accepted group %d %d\n", this->scope->name, this->rank, group->i,group->j);
        recvlist.push_back(group);
        group->owner = ambient::rank();
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
    void assign(void_spt* ptr, int i, int j, int k)
    {
        asmp.assign(ptr->group(i, j, k));
    }

}
