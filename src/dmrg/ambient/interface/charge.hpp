#ifndef AMBIENT_INTERFACE_CHARGE_H
#define AMBIENT_INTERFACE_CHARGE_H

namespace ambient {

    class charge { // workload of individual rank in terms of workgroups 
    public:
        charge():accept(false){ this->set_smp_group("ambient"); }
        charge& operator()(const int rank)
        {
            if(rank == ambient::rank(this->smp_group)){ accept = true; }
            else{ accept = false; target = rank; }
            return *this;
        }
        charge& operator+=(workgroup* group)
        {
            if(accept){
                accept = false;
                printf("R%d: I've accepted group %d %d\n", ambient::rank(this->smp_group), group->i,group->j);
                recvlist.push_back(group);
                group->owner = ambient::rank(this->smp_group);
            }else if(group->owner == ambient::rank(this->smp_group)){
                group->owner = target;
                sendlist.push_back(group);
            }
        }
        void set_smp_group(const char* smp_group)
        {
            this->smp_group = smp_group;
        }
        const char* smp_group;
        bool accept;
        int target;
        std::list<workgroup*> sendlist;
        std::list<workgroup*> recvlist;
    } charge;
}
#endif
