#ifndef AMBIENT_CHANNELS_PACKETS
#define AMBIENT_CHANNELS_PACKETS
#include <mpi.h>
#include "ambient/channels/ichannel.h"
#include "ambient/channels/packets/packet_t.h"

namespace ambient { namespace channels {

    class packet : public ichannel::packet
    {
    public:
        void* data;
        MPI_Datatype mpi_t;
        int lifetime;
        const packet_t& type;
        const packet_t& get_t();
        MPI_Datatype get_mpi_t();
        int   get_t_code();
        const void* get(int field);
        size_t get_bound(size_t field);
        void* get_memory();
        bool  disposable();

        template<typename T>
        T get(int field){
            return *(T*)this->get(field);
        }

        void set(int field, const void* value);
        void set(int field, int value);
        packet(const packet_t& type, const void* memory);
        packet(const packet_t& type, void* memory, va_list& fields); // used in auxiliary.hpp
        packet(const packet& p);
       ~packet();
    };

    packet* pack(const packet_t& type, void* memory, ...);
    packet* unpack(const packet_t& type, void* memory);

} }

#endif
