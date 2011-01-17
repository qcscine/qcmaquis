#ifndef AMBIENT_PACKETS_PACKET_TYPE_H
#define AMBIENT_PACKETS_PACKET_TYPE_H
#include <memory.h>
#include <mpi.h>

#define LEN(x) sizeof(x)/sizeof(int)
#define PACK   static int sizes[] = 
#define FIELDS friend class packet_t; protected: MPI_Datatype
#define BASE_FIELDS     MPI_Datatype
#define __A_PACKET__    MPI_Datatype& __a_packet_start = 
#define __A_CODE(code)  construct(code, LEN(sizes), \
                        sizes, &__a_packet_start);
namespace ambient{ namespace packets{

    class packet_t
    {
// PACKET-TYPE WORK LOGIC //////////////////////////////////////////////////////////////
    private:
        packet_t(packet_t const&){};             // copy constructor is private
        packet_t& operator=(packet_t const&){};  // assignment operator is private
    public:
        template<typename T>
        static T& get(){
            static T* singleton = NULL;
            if(!singleton){ 
                singleton = new T(); 
                type_map(singleton->t_code, (packet_t*)singleton);
            }
            return *singleton;
        }
        static packet_t& type_map(char t_code, const packet_t* type = NULL){
            static std::map<char,const packet_t*> map;
            if(type != NULL) map.insert(std::pair<char,const packet_t*>(t_code,type));
            else return const_cast<packet_t&>(*(map.find(t_code)->second));
        }
        void change_field_size(int field, int size);
        void fill_packet(void* memory, char type, va_list& fields) const;
        void commit();
    protected:
        void construct(char code, int count, const int* sizes, const MPI_Datatype* types);
    public:
        char t_code;
        int  t_size;                                   // memory size of the type
        int  count;                                    // number of types used
        int* sizes;                                    // block-lengths of used types
        MPI_Aint* displacements;                       // memory placements of type-blocks
        MPI_Datatype  mpi_t;                           // resulting MPI Datatype
        MPI_Datatype* compounds;                       // MPI compound datatypes
        BASE_FIELDS type;                              // mandatory packet-type first field
    protected:
        packet_t();
    };

} }
#endif
