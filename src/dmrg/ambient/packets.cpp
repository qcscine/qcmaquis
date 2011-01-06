#include "ambient/packets.h"

void packet_type::type_displacements(int count, MPI_Datatype* types, int* sizes)
{
    MPI_Aint extent;
    this->displacements = (MPI_Aint*)realloc(this->displacements, (this->count+count)*sizeof(MPI_Aint));
    this->sizes = (int*)realloc(this->sizes, (this->count+count)*sizeof(int));
    memcpy(&(this->sizes[this->count]), sizes, count*sizeof(int));

    for(int i = 0; i < count; i++)
    {
        this->displacements[this->count + i] = this->size;
        MPI_Type_extent(types[i], &extent);
        this->size += extent*sizes[i];
    }
    this->count += count;
}

void packet_type::change_field_size(int field, int size)
{
    MPI_Aint extent;
    MPI_Type_extent((&this->type)[field], &extent);
    int delta = extent*(size - this->sizes[field]);
    this->sizes[field] = size;

    for(int i = field+1; i < this->count; i++)
        this->displacements[i] += delta;
}

// note: everything except for INT is passed by pointers
void packet_type::instance(void* memory, ...)
{
    void* iterator;
    void* source;
    int   type_size;
    va_list fields;
    MPI_Datatype* types = &this->type;
    va_start(fields, memory); 
    for(int i = 0; i < this->count; i++)
    {
        iterator = (void*)((size_t)memory + this->displacements[i]);
        if(sizes[i] == 1 && types[i] == MPI_INT)
            *((int *)iterator) = va_arg(fields, int);
        else if((source = va_arg(fields, void*)) != NULL){
            MPI_Type_size(types[i], &type_size);
            memcpy(iterator, source, sizes[i]*type_size);
        }
    }
    va_end(fields);
}
