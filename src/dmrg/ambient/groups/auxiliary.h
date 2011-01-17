#ifndef AMBIENT_GROUPS_AUX_H
#define AMBIENT_GROUPS_AUX_H

namespace ambient{ namespace groups{

    template<typename T>
    packet* recv(void* memory){
        packet_manager::instance().recv(get_t<T>(), memory);
        return unpack<T>(memory);
    }

    void send(packet* pack, int dest = -1)
    {
        if(dest != -1){
            packet_manager::instance().send(pack, dest);
        }else{
            if(pack->get_t().compounds[1] != MPI_INT) printf("Warning: the dest field (#1) is not of type MPI_INT!\n");
            packet_manager::instance().send(pack, *(int*)pack->get(A_DEST_FIELD));
        }
    }

} }
#endif
