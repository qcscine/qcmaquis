#include "ambient/groups/packet_manager.h"
#include "ambient/groups/group.h"

namespace ambient{ namespace groups{

    packet_manager::packet_manager(MPI_Comm* comm){
        this->comm = comm;
    };

    void packet_manager::send(packet* pack, int dest)
    {
        if(pack->get_mpi_t() != pack->mpi_t) printf("Erroneous packet type: the type has changed since packet's creation\n");
        MPI_Send(pack->data, 1, pack->mpi_t, dest, pack->get_t_code(), *this->comm);
    }

    void packet_manager::recv(const packet_t& type, void* memory)
    {
        MPI_Recv(memory, 1, type.mpi_t, MPI_ANY_SOURCE, type.t_code, *this->comm, MPI_STATUS_IGNORE);
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

    void send(packet* pack, const char* grp, int dest = -1)
    {
        send(pack, group_map(grp), dest);
    }

    packet* recv(const packet_t& type, group* grp, void* memory){
        grp->manager->recv(type, memory);
        return unpack(type, memory);
    }

    packet* recv(const packet_t& type, const char* grp, void* memory){
        return recv(type, group_map(grp), memory);
    }

} }
