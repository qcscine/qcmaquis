#ifndef AMBIENT_PACKETS_H
#define AMBIENT_PACKETS_H
#include <mpi.h>
#include "ambient/packets/packet_t.h"

namespace ambient{ namespace packets{

    class packet
    {
    public:
        void* data;
        MPI_Datatype mpi_t;

        packet_t* get_t();
        MPI_Datatype get_mpi_t();
        char  get_t_code();
        void* get(int field);
        void  set(int field, void* value);
        void  set(int field, int value);
        void  send(int dest = -1);
        packet(packet_t* type, void* memory, ...);
        packet(packet_t* type, ...);
        packet(void* memory);

        packet(packet_t* type, void* memory, va_list& fields);
    };

} }
#endif
