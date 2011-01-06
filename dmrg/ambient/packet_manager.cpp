#include "ambient/packet_manager.h"

packet_manager* packet_manager::instance(){
    static packet_manager* singleton = NULL;
    if(!singleton) singleton = new packet_manager();
    return singleton;
}
packet_manager::packet_manager(){};

void packet_manager::commit(packet_type* type)
{
    MPI_Type_struct(type->count,                            // count - number of types
                    type->sizes,                            // number of instances of the same type (blocklengths)
                    type->displacements,                    // starting positions
                    &type->type,                            // types being used
                    &type->mpi_type);

    MPI_Type_commit(&type->mpi_type);
}
