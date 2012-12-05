#include "ambient/utils/timings.hpp"

namespace ambient { namespace channels { namespace mpi {

    inline packet_t::packet_t()
    : t_size(0), count(0), displacements(NULL), sizes(NULL), compounds(NULL), mpi_t(MPI_DATATYPE_NULL)
    {
        type = MPI_INT;
        int sizes[] = { 1 };
        construct('0', 1, sizes, &type);
    }

    inline void packet_t::construct(int code, int count, const int* sizes, const MPI_Datatype* types){
        MPI_Aint extent;
        this->t_code = code;
        this->compounds = (MPI_Datatype*)realloc(this->compounds, (this->count+count)*sizeof(MPI_Datatype));
        this->displacements = (MPI_Aint*)realloc(this->displacements, (this->count+count)*sizeof(MPI_Aint));
        this->sizes = (int*)realloc(this->sizes, (this->count+count)*sizeof(int));
        memcpy(&(this->sizes[this->count]), sizes, count*sizeof(int));
        memcpy(&(this->compounds[this->count]), types, count*sizeof(MPI_Datatype));

        for(int i = 0; i < count; i++)
        {
            this->displacements[this->count + i] = this->t_size;
            MPI_Type_extent(types[i], &extent);
            this->t_size += extent*sizes[i];
        }
        this->count += count;
    }

    inline void packet_t::change_field_size(int field, int size){
        MPI_Aint extent;
        MPI_Type_extent(this->compounds[field], &extent);
        int delta = extent*(size - this->sizes[field]);
        this->sizes[field] = size;

        for(int i = field+1; i < this->count; i++)
            this->displacements[i] += delta;
        this->t_size += delta;
    }

    // note: everything except for INT is passed by pointers
    inline void packet_t::fill_packet(void* memory, int type, va_list& fields) const {
        void* iterator;
        void* source;
        int   t_size;
        *(int*)((size_t)memory + this->displacements[0]) = type; // type
        //*(int*)((size_t)memory + this->displacements[1]) = 0;  // see alloc_t

        for(int i = 2; i < this->count; i++)
        {
            iterator = (void*)((size_t)memory + this->displacements[i]);
            if(sizes[i] == 1 && this->compounds[i] == MPI_INT)
                *((int *)iterator) = va_arg(fields, int);
            else if((source = va_arg(fields, void*)) != NULL){
                MPI_Type_size(this->compounds[i], &t_size);
                memcpy(iterator, source, sizes[i]*t_size);
            }
        }
    }

    inline size_t packet_t::get_size() const {
        return this->t_size;
    }

    inline size_t packet_t::get_bound(size_t field) const {
        return this->displacements[field];
    }

    inline void packet_t::commit(){
        MPI_Type_create_struct(this->count,            // count - number of types
                               this->sizes,            // number of instances of the same type (blocklengths)
                               this->displacements,    // starting positions
                               this->compounds,        // types being used
                               &this->mpi_t);

        MPI_Type_commit(&this->mpi_t);
    }

    // {{{ memory management (C-style intrusive ptr)

    inline void* alloc_t(const packet_t& type){
        void* memory = malloc(type.get_size());
        memset(memory, 0, type.get_size());
        return memory;
    }

    // }}}
    
} } }
