#include "ambient/groups/packet_manager.h"

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

} }
