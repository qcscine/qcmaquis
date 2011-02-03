#ifndef AMBIENT_INTERFACE_CHARGE_H
#define AMBIENT_INTERFACE_CHARGE_H

#include "ambient/interface/workgroup.h"

namespace ambient {

    class charge { // workload of individual rank in terms of workgroups 
    public:
        charge();
        charge& operator()(const int rank);
        charge& operator+=(workgroup* group);
        void set_smp_group(const char* smp_group);
        const char* smp_group;
        bool accept;
        int target;
        std::list<workgroup*> sendlist;
        std::list<workgroup*> recvlist;
    };

    extern charge workload;
}
#endif
