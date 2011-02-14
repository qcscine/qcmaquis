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
        if(asmp.rank == UNDEFINED_RANK) return;
        workgroup* group = ptr->group(i,j,k);
//        printf("%s: p%d: I've accepted group %d %d of id%d\n", asmp.get_scope()->name, asmp.rank, group->i, group->j, (*(group->profile))->id );
        group->owner = ambient::rank();
        if(!group->initialized){
            group->data = ptr->init_fp(group);
            group->initialized = true;
        }
        ptr->layout->add_segment_entry(ambient::rank(), i, j, k);
    }

}
