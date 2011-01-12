#include "ambient/packet_manager.h"

namespace ambient
{
    packet_manager* packet_manager::instance(){
        static packet_manager* singleton = NULL;
        if(!singleton) singleton = new packet_manager();
        return singleton;
    }
    packet_manager::packet_manager(){};

    void packet_manager::send(packet* pack, int dest)
    {
// correctness check:
        if(pack->get_mpi_t() != pack->mpi_type) printf("Erroneous packet type: the type has changed since packet's creation\n");
        MPI_Send(pack->data, 1, pack->mpi_type, dest, pack->get_t_code(), MPI_COMM_WORLD);
    }
}
