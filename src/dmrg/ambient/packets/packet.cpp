#include "ambient/packets/types.h"
#include "ambient/packets/packet.h"
#include <stdlib.h>
#include <stdarg.h>

namespace ambient{ namespace packets {

    packet::packet(const packet_t& type, void* memory, ...){
        va_list fields;
        va_start(fields, memory); 
        type.fill_packet(memory, type.t_code, fields);
        va_end(fields);
        this->data = memory;
        this->mpi_t = this->get_mpi_t(); // handy errors checking
    }

    packet::packet(const packet_t& type, ...){
        va_list fields;
        void* memory;
        va_start(fields, type); 
        memory = malloc(type.t_size);
        type.fill_packet(memory, type.t_code, fields);
        va_end(fields);
        this->data = memory;
        this->mpi_t = this->get_mpi_t(); // handy errors checking
    }

    packet::packet(const packet_t& type, void* memory, va_list& fields){
        type.fill_packet(memory, type.t_code, fields);
        this->data = memory;
        this->mpi_t = this->get_mpi_t(); // handy errors checking
    }
 
    packet::packet(const void* memory){
        this->data = const_cast<void*>(memory);
        this->mpi_t = this->get_mpi_t(); // handy errors checking
    }

    packet_t& packet::get_t()
    {
        return packet_t::type_map(this->get_t_code());
    }

    char packet::get_t_code()
    {
        return *(char*)this->data;
    }

    MPI_Datatype packet::get_mpi_t()
    {
        return this->get_t().mpi_t;
    }

    const void* packet::get(int field)
    {
        return (void*)((size_t)this->data + this->get_t().displacements[field]);
    }

    void packet::set(int field, const void* value)
    {
        int   t_size;
        if(field == A_TYPE_FIELD) printf("Warning: attempting to modify the type of the packet.\n");
        void* ptr = const_cast<void*>(this->get(field));
        MPI_Type_size(this->get_t().compounds[field], &t_size);
        memcpy(ptr, value, this->get_t().sizes[field]*t_size);
    }

    void packet::set(int field, int value)
    {
        if(field == A_TYPE_FIELD) printf("Warning: attempting to modify the type of the packet.\n");
        int* ptr = (int*)const_cast<void*>(this->get(field));
        *ptr = value;
    }

    packet* pack(const packet_t& type, void* memory, ...){
        packet* instance;
        va_list fields;
        va_start(fields, memory); 
        instance = new packet(type, memory, fields);
        va_end(fields);
        return instance;
    }

    packet* unpack(const packet_t& type, void* memory){
        *(char*)memory = type.t_code;
        return new packet(memory);
    }

} }
