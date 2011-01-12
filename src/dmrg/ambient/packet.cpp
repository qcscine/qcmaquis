#include "ambient/types.h"
#include "ambient/packet.h"
#include "ambient/packet_manager.h"

packet::packet(packet_type* type, void* memory, ...){
    va_list fields;
    va_start(fields, memory); 
    type->fill_packet(memory, type->type_code, fields);
    va_end(fields);
    this->data = memory;
    this->mpi_type = this->get_mpi_type(); // handy errors checking
}

packet::packet(packet_type* type, ...){
    va_list fields;
    void* memory;
    va_start(fields, type); 
    memory = malloc(type->type_size);
    type->fill_packet(memory, type->type_code, fields);
    va_end(fields);
    this->data = memory;
    this->mpi_type = this->get_mpi_type(); // handy errors checking
}
    
packet::packet(void* memory){
    this->data = memory;
    this->mpi_type = this->get_mpi_type(); // handy errors checking
}

packet_type* packet::get_type()
{
    return packet_type::type_map(*(char*)this->data);
}

char packet::get_type_code()
{
    return this->get_type()->type_code;
}

MPI_Datatype packet::get_mpi_type()
{
    return this->get_type()->mpi_type;
}

void* packet::get(int field)
{
    return (void*)((size_t)this->data + this->get_type()->displacements[field]);
}

void packet::set(int field, void* value)
{
    int   type_size;
    if(field == A_TYPE_FIELD) printf("Warning: attempting to modify the type of the packet.\n");
    void* ptr = this->get(field);
    MPI_Type_size(this->get_type()->compounds[field], &type_size);
    memcpy(ptr, value, this->get_type()->sizes[field]*type_size);
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
        if(this->get_type()->compounds[1] != MPI_INT) printf("Warning: the dest field (#1) is not of type MPI_INT!\n");
        packet_manager::instance()->send(this, *(int*)this->get(A_DEST_FIELD));
    }
}
