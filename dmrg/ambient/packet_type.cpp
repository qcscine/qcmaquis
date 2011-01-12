#include "ambient/packet_type.h"

packet_type::packet_type(): type_size(0), count(0), 
                            displacements(NULL), sizes(NULL),
                            compounds(NULL)
{
    type = MPI_BYTE;
    int sizes[] = { 1 };
    construct('0', 1, sizes, &type);
}

void packet_type::construct(char code, int count, int* sizes, MPI_Datatype* types)
{
    MPI_Aint extent;
    this->type_code = code;
    this->compounds = (MPI_Datatype*)realloc(this->compounds, (this->count+count)*sizeof(MPI_Datatype));
    this->displacements = (MPI_Aint*)realloc(this->displacements, (this->count+count)*sizeof(MPI_Aint));
    this->sizes = (int*)realloc(this->sizes, (this->count+count)*sizeof(int));
    memcpy(&(this->sizes[this->count]), sizes, count*sizeof(int));
    memcpy(&(this->compounds[this->count]), types, count*sizeof(MPI_Datatype));

    for(int i = 0; i < count; i++)
    {
        this->displacements[this->count + i] = this->type_size;
        MPI_Type_extent(types[i], &extent);
        this->type_size += extent*sizes[i];
    }
    this->count += count;
}

void packet_type::change_field_size(int field, int size)
{
    MPI_Aint extent;
    MPI_Type_extent(this->compounds[field], &extent);
    int delta = extent*(size - this->sizes[field]);
    this->sizes[field] = size;

    for(int i = field+1; i < this->count; i++)
        this->displacements[i] += delta;
    this->type_size += delta;
}

// note: everything except for INT is passed by pointers
void packet_type::fill_packet(void* memory, const char type, va_list& fields)
{
    void* iterator;
    void* source;
    int   type_size;
    *(char*)memory = type;

    for(int i = 1; i < this->count; i++)
    {
        iterator = (void*)((size_t)memory + this->displacements[i]);
        if(sizes[i] == 1 && this->compounds[i] == MPI_INT)
            *((int *)iterator) = va_arg(fields, int);
        else if((source = va_arg(fields, void*)) != NULL){
            MPI_Type_size(this->compounds[i], &type_size);
            memcpy(iterator, source, sizes[i]*type_size);
        }
    }
}

void packet_type::commit()
{
    MPI_Type_create_struct(this->count,                            // count - number of types
                           this->sizes,                            // number of instances of the same type (blocklengths)
                           this->displacements,                    // starting positions
                           this->compounds,                        // types being used
                           &this->mpi_type);

    MPI_Type_commit(&this->mpi_type);
}
