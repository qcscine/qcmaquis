#include "ambient/channels/packets/types.h"
#include "ambient/channels/packets/packet.h"
#include <stdlib.h>
#include <stdarg.h>

namespace ambient { namespace channels {

    packet::packet(const packet_t& type, void* memory, va_list& fields) 
    : type(type), lifetime(1)
    {
        type.fill_packet(memory, type.t_code, fields);
        this->data = memory;
        this->mpi_t = this->get_mpi_t();
    }
 
    packet::packet(const packet_t& type, const void* memory)
    : type(type), lifetime(1)
    {
        this->data = const_cast<void*>(memory);
        this->mpi_t = this->get_mpi_t();
    }

    packet::packet(const packet& p) 
    : type(p.type)
    {
        memcpy(this->data, p.data, p.type.get_size());
        reset_usage(this->data);
        this->mpi_t = p.mpi_t;
    }

    packet::~packet(){
    }

    bool packet::disposable(){
        return (this->lifetime > 0 && (--this->lifetime) == 0);
    }

    const packet_t& packet::get_t(){
        return this->type;
    }

    MPI_Datatype packet::get_mpi_t(){
        return this->type.mpi_t;
    }

    int packet::get_t_code(){
        return *(int*)this->data;
    }

    const void* packet::get(int field){
        return (void*)((size_t)this->data + this->get_t().displacements[field]);
    }

    size_t packet::get_bound(size_t field){
        return this->get_t().displacements[field];
    }

    void* packet::get_memory(){
        return this->data;
    }

    void packet::set(int field, const void* value){
        int   t_size;
        if(field == A_TYPE_FIELD) printf("Warning: attempting to modify the type of the packet.\n");
        void* ptr = const_cast<void*>(this->get(field));
        MPI_Type_size(this->get_t().compounds[field], &t_size);
        memcpy(ptr, value, this->get_t().sizes[field]*t_size);
    }

    void packet::set(int field, int value){
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

    packet* unpack(const packet_t& type, void* memory){ // stack (dies quickly)
        *(int*)memory = type.t_code;
        return new packet(type, memory);
    }

} }
