#ifndef AMBIENT_CHANNELS_MPI_PACKETS_AUX
#define AMBIENT_CHANNELS_MPI_PACKETS_AUX

namespace ambient { namespace channels { namespace mpi {

    template<typename T>
    T& get_t(){
        return packet_t::get<T>();
    }

    template<typename T>
    MPI_Datatype get_mpi_t(){
        return get_t<T>().mpi_t;
    }

    template<typename T>
    size_t sizeof_t(int field = -1){
        if(field == -1) return get_t<T>().get_size();
        else return get_t<T>().sizes[field];
    }

    template<typename T>
    void commit_t(){
        if(get_mpi_t<T>() != MPI_DATATYPE_NULL) MPI_Type_free(&get_t<T>().mpi_t);
        get_t<T>().commit();
    }

    template<typename T>
    void change_t(int field, int size){
        get_t<T>().change_field_size(field, size);
        commit_t<T>();
    }

    template<typename T>
    void* alloc_t(){
        return alloc_t(get_t<T>());
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

    template<typename T>
    packet* unpack(void* memory){
        return unpack(get_t<T>(), memory);
    }

} } }

#endif
