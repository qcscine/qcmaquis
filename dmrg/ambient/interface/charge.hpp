#ifndef AMBIENT_INTERFACE_CHARGE_H
#define AMBIENT_INTERFACE_CHARGE_H

namespace ambient {

    class charge { // workload of individual rank in terms of workgroups 
    public:
        charge():accept(false){ }
        charge& operator()(const int rank)
        {
            if(rank == 0){ accept = true; }
            else{ accept = false; target = rank; }
            return *this;
        }
        charge& operator+=(workgroup* group)
        {
            if(accept){
                accept = false;
                zout << "I've accepted group " << group->i << " " << group->j << std::endl;
                recvlist.push_back(group);
                group->owner = 0;
            }else if(group->owner == 0){
                group->owner = target;
                sendlist.push_back(group);
            }
        }
        bool accept;
        int target;
        std::list<workgroup*> sendlist;
        std::list<workgroup*> recvlist;
    } charge;
}
#endif
