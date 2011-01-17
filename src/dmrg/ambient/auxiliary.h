#ifndef AMBIENT_AUX_H
#define AMBIENT_AUX_H

namespace ambient{ namespace groups{

    template<typename T>
    packet* recv(group* grp, void* memory){
        grp->manager->recv(get_t<T>(), memory);
        return unpack<T>(memory);
    }

    void send(packet* pack, group* grp, int dest = -1)
    {
        if(dest != -1){
            grp->manager->send(pack, dest);
        }else{
            if(pack->get_t().compounds[1] != MPI_INT) printf("Warning: the dest field (#1) is not of type MPI_INT!\n");
            grp->manager->send(pack, *(int*)pack->get(A_DEST_FIELD));
        }
    }

} }
#endif
