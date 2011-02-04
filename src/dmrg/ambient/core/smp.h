#ifndef AMBIENT_INTERFACE_CHARGE_H
#define AMBIENT_INTERFACE_CHARGE_H

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
        void set_smp_group(const char* smp_group);
        const char* smp_group;
        std::list<workgroup*> sendlist;
        std::list<workgroup*> recvlist;
        int id;
    };
    void assign(p_profile_s* ptr, int i, int j = 0, int k = 0);

}
#endif
