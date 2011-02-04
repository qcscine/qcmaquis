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
    smp::smp(){ }

    smp& smp::operator()(const int rank)
    {
        return *this;
    }
    void smp::assign(workgroup* group)
    {
        printf("%s: p%d: I've accepted group %d %d\n", this->smp_group, this->id, group->i,group->j);
        recvlist.push_back(group);
        group->owner = ambient::rank(this->smp_group);
    }
    void smp::set_smp_group(const char* smp_group)
    {
        this->smp_group = smp_group;
        this->id = ambient::rank(smp_group);
    }

    void assign(void_spt* ptr, int i, int j, int k)
    {
        if(asmp.id == UNDEFINED_ID) return;
        asmp.assign(ptr->group(i, j, k));
    }

}
