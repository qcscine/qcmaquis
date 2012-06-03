#ifndef AMBIENT_CHANNELS_MPI_PACKETS_INSTANCE
#define AMBIENT_CHANNELS_MPI_PACKETS_INSTANCE

namespace ambient { namespace channels { namespace mpi {

    class packet
    {
    public:
        void* data;
        MPI_Datatype mpi_t;
        int lifetime;
        const packet_t& type;
        inline const packet_t& get_t();
        inline MPI_Datatype get_mpi_t();
        inline int   get_t_code();
        inline const void* get(int field);
        inline size_t get_bound(size_t field);
        inline void* get_memory();
        inline bool  disposable();

        template<typename T>
        inline T get(int field){
            return *(T*)this->get(field);
        }

        inline void set(int field, const void* value);
        inline void set(int field, int value);
        inline packet(const packet_t& type, const void* memory);
        inline packet(const packet_t& type, void* memory, va_list& fields); // used in auxiliary.hpp
        inline packet(const packet& p);
        inline ~packet();
    };

    inline packet* pack(const packet_t& type, void* memory, ...);
    inline packet* unpack(const packet_t& type, void* memory);

} } }

#endif
