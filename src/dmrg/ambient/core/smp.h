#ifndef AMBIENT_CORE_SMP_H
#define AMBIENT_CORE_SMP_H

#include "ambient/core/operation.h"
#include "ambient/core/p_profile.h"
#include "ambient/core/workgroup.h"
#include "ambient/core/layout.h"

namespace ambient {

typedef ambient::p_profile   void_pt;
typedef ambient::p_profile_s void_spt;

    class smp { // workload of individual rank in terms of workgroups 
    private:
        smp();
        smp(smp const&);
        smp& operator=(smp const&);
    public:
        static smp& instance();

    public:
        smp& operator()(const int rank);
        void set_scope(groups::group* scope);
        void set_scope(const char* scope);
        groups::group* get_scope();
        int scope_size;
        int rank;
        core::operation* op;
    private:
        groups::group* scope;
        core::layout_table* assignment;
        std::list<workgroup*>  sendlist;
        std::list<workgroup*>  recvlist;
    };
    void assign(void_spt* ptr, int i, int j = 0, int k = 0);

}
#endif
