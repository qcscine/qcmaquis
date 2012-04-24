#ifndef AMBIENT_CHANNELS_PACKETS_BASETYPE
#define AMBIENT_CHANNELS_PACKETS_BASETYPE
#include <memory.h>
#include <mpi.h>
#include <stdlib.h>
#include <map>

#define LEN(x) sizeof(x)/sizeof(int)
#define BASE_FIELDS       MPI_Datatype
#define __a_pack          int a_packet_types_sizes[] = 
#define __a_fields__      friend class packet_t; protected: MPI_Datatype
#define __a_flex_fields__ friend class packet_t; public: MPI_Datatype
#define __a_packet__      MPI_Datatype& __a_packet_start = 
#define __a_code(code)    construct(code, LEN(a_packet_types_sizes), \
                          a_packet_types_sizes, &__a_packet_start);
namespace ambient { namespace channels {

    class packet_t
    {
// PACKET-TYPE WORK LOGIC //////////////////////////////////////////////////////////////
    private:
        packet_t(packet_t const&);             // copy constructor is private
        packet_t& operator=(packet_t const&);  // assignment operator is private
    public:
        template<typename T>
        static T& get(){
            static T* singleton = NULL;
            if(!singleton){ 
                singleton = new T();
                singleton->commit(); 
                type_map(singleton->t_code, (packet_t*)singleton);
            }
            return *singleton;
        }
        static packet_t& type_map(int t_code, const packet_t* type = NULL){
            static std::map<int,const packet_t*> map;
            if(type != NULL) map.insert(std::pair<int,const packet_t*>(t_code,type));
            return const_cast<packet_t&>(*(map.find(t_code)->second));
        }
        void change_field_size(int field, int size);
        void fill_packet(void* memory, int type, va_list& fields) const;
        size_t get_size() const;
        size_t get_bound(size_t field) const;
        void commit();
    protected:
        void construct(int code, int count, const int* sizes, const MPI_Datatype* types);
    public:
        int  t_code;
        int  t_size;                                   // memory size of the type
        int  count;                                    // number of types used
        int* sizes;                                    // block-lengths of used types
        MPI_Aint* displacements;                       // memory placements of type-blocks
        MPI_Datatype  mpi_t;                           // resulting MPI Datatype
        MPI_Datatype* compounds;                       // MPI compound datatypes
        BASE_FIELDS type, usage;                       // mandatory packet-type first field
    protected:
        packet_t();
    };

    void* alloc_t(const packet_t& type);
    void reset_usage(void* memory);
    void checked_free(void* memory);
    void use(void* memory);
    void unuse(void* memory);

} }

#endif
