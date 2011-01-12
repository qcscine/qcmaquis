#ifndef AMBIENT_AUXILIARY_H
#define AMBIENT_AUXILIARY_H

namespace ambient
{
    template<typename T>
    T* get_t(){
        return packet_t::get<T>();
    }

    template<typename T>
    MPI_Datatype get_mpi_t(){
        return get_t<T>()->mpi_t;
    }

    template<typename T>
    size_t sizeof_t(int field = -1){
        if(field == -1) return get_t<T>()->t_size;
        else return get_t<T>()->sizes[field];
    }

    template<typename T>
    void change_t(int field, int size){
        get_t<T>()->change_field_size(field, size);
    }

    template<typename T>
    void commit_t(){
        get_t<T>()->commit();
    }

    template<typename T>
    packet* pack(void* memory, ...){
        packet* instance;
        va_list fields;
        va_start(fields, memory); 
        instance = new packet(get_t<T>(), memory, fields);
        va_end(fields);
        return instance;
    }

    packet* unpack(void* memory){ 
        return new packet(memory);
    }

    template<typename T>
    packet* unpack(void* memory){
        *(char*)memory = get_t<T>()->t_code;
        return unpack(memory);
    }
}

#endif
