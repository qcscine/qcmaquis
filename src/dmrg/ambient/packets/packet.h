#ifndef AMBIENT_PACKETS_PACKET_H
#define AMBIENT_PACKETS_PACKET_H
#include <mpi.h>
#include "ambient/packets/packet_t.h"

namespace ambient{ namespace packets{

    class packet
    {
    public:
        void* data;
        MPI_Datatype mpi_t;

        const packet_t& type;
        const packet_t& get_t();
        MPI_Datatype get_mpi_t();
        char  get_t_code();
        const void* get(int field);

        template<typename T>
        T get(int field){
            return *(T*)this->get(field);
        }

        void set(int field, const void* value);
        void set(int field, int value);
        packet(const packet_t& type, const void* memory);
        packet(const packet_t& type, void* memory, va_list& fields); // used in auxiliary.hpp
    };

    packet* pack(const packet_t& type, void* memory, ...);
    packet* unpack(const packet_t& type, void* memory);

} }
#endif
