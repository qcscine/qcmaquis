#include "ambient/groups/packet_manager.h"

namespace ambient{ namespace groups{

    packet_manager& packet_manager::instance()
    {
        static packet_manager* singleton = NULL;
        if(!singleton) singleton = new packet_manager();
        return *singleton;
    }
    packet_manager::packet_manager(){};

    void packet_manager::set_comm(MPI_Comm comm)
    {
        this->comm = comm;
    }

    void packet_manager::send(packet* pack, int dest)
    {
// correctness check:
        if(pack->get_mpi_t() != pack->mpi_t) printf("Erroneous packet type: the type has changed since packet's creation\n");
        MPI_Send(pack->data, 1, pack->mpi_t, dest, pack->get_t_code(), this->comm);
    }

    void packet_manager::recv(const packet_t& type, void* memory)
    {
        MPI_Recv(memory, 1, type.mpi_t, MPI_ANY_SOURCE, type.t_code, this->comm, MPI_STATUS_IGNORE);
    }

} }
