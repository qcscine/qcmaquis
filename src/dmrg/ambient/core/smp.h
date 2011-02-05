#ifndef AMBIENT_INTERFACE_CHARGE_H
#define AMBIENT_INTERFACE_CHARGE_H

#include "ambient/core/operation.h"
#include "ambient/interface/p_profile.h"
#include "ambient/interface/workgroup.h"

namespace ambient {

    class smp { // workload of individual rank in terms of workgroups 
    private:
        smp();
        smp(smp const&);
        smp& operator=(smp const&);
    public:
        static smp& instance();

    public:
        smp& operator()(const int rank);
        void assign(workgroup* group);
        void set_scope(const groups::group* scope);
        void set_scope(const char* scope);
        void get_info(core::operation* op);
        const groups::group* scope;
        int scope_size;
        int id;
        bool accept;
    private:
        std::list<workgroup*> sendlist;
        std::list<workgroup*> recvlist;
    };
    void assign(p_profile_s* ptr, int i, int j = 0, int k = 0);

}
#endif
