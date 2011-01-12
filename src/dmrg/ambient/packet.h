#ifndef AMBIENT_PACKETS_H
#define AMBIENT_PACKETS_H
#include <mpi.h>
#include "ambient/packet_type.h"

    
class packet
{
public:
    void* data;
    MPI_Datatype mpi_type;

    packet_type* get_type();
    MPI_Datatype get_mpi_type();
    char  get_type_code();
    void* get(int field);
    void  set(int field, void* value);
    void  set(int field, int value);
    void  send(int dest = -1);
    packet(packet_type* type, void* memory, ...);
    packet(packet_type* type, ...);
    packet(void* memory);
};

#endif
