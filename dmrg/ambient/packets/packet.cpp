#include "ambient/packets/types.h"
#include "ambient/packets/packet.h"
#include "ambient/packets/packet_manager.h"

namespace ambient{ namespace packets {

    packet::packet(packet_t* type, void* memory, ...){
        va_list fields;
        va_start(fields, memory); 
        type->fill_packet(memory, type->t_code, fields);
        va_end(fields);
        this->data = memory;
        this->mpi_t = this->get_mpi_t(); // handy errors checking
    }

    packet::packet(packet_t* type, ...){
        va_list fields;
        void* memory;
        va_start(fields, type); 
        memory = malloc(type->t_size);
        type->fill_packet(memory, type->t_code, fields);
        va_end(fields);
        this->data = memory;
        this->mpi_t = this->get_mpi_t(); // handy errors checking
    }

    packet::packet(packet_t* type, void* memory, va_list& fields){
        type->fill_packet(memory, type->t_code, fields);
        this->data = memory;
        this->mpi_t = this->get_mpi_t(); // handy errors checking
    }
 
    packet::packet(void* memory){
        this->data = memory;
        this->mpi_t = this->get_mpi_t(); // handy errors checking
    }

    packet_t* packet::get_t()
    {
        return packet_t::type_map(*(char*)this->data);
    }

    char packet::get_t_code()
    {
        return this->get_t()->t_code;
    }

    MPI_Datatype packet::get_mpi_t()
    {
        return this->get_t()->mpi_t;
    }

    void* packet::get(int field)
    {
        return (void*)((size_t)this->data + this->get_t()->displacements[field]);
    }

    void packet::set(int field, void* value)
    {
        int   t_size;
        if(field == A_TYPE_FIELD) printf("Warning: attempting to modify the type of the packet.\n");
        void* ptr = this->get(field);
        MPI_Type_size(this->get_t()->compounds[field], &t_size);
        memcpy(ptr, value, this->get_t()->sizes[field]*t_size);
    }

    void packet::set(int field, int value)
    {
        if(field == A_TYPE_FIELD) printf("Warning: attempting to modify the type of the packet.\n");
        int* ptr = (int*)this->get(field);
        *ptr = value;
    }

    void packet::send(int dest)
    {
        if(dest != -1){
            packet_manager::instance()->send(this, dest);
        }else{
            if(this->get_t()->compounds[1] != MPI_INT) printf("Warning: the dest field (#1) is not of type MPI_INT!\n");
            packet_manager::instance()->send(this, *(int*)this->get(A_DEST_FIELD));
        }
    }

} }
