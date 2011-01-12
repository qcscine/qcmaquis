#ifndef AMBIENT_PACKET_TYPE_H
#define AMBIENT_PACKET_TYPE_H
#include <memory.h>
#include <mpi.h>

#define LEN(x) sizeof(x)/sizeof(int)
#define PACK   static int sizes[] = 
#define FIELDS friend class packet_type; protected: MPI_Datatype
#define BASE_FIELDS     MPI_Datatype
#define __A_PACKET__    MPI_Datatype& __a_packet_start = 
#define __A_CODE(code)  construct(code, LEN(sizes), \
                        sizes, &__a_packet_start);

class packet_type
{
// PACKET-TYPE WORK LOGIC //////////////////////////////////////////////////////////////
private:
    packet_type(packet_type const&){};             // copy constructor is private
    packet_type& operator=(packet_type const&){};  // assignment operator is private
public:
    template<class typeT>
    static typeT* get(){
        static typeT* singleton = NULL;
        if(!singleton){ 
            singleton = new typeT(); 
            type_map(singleton->type_code, (packet_type*)singleton);
        }
        return singleton;
    }
    static packet_type* type_map(char type_code, packet_type* type = NULL){
        static std::map<char,packet_type*> map;
        if(type != NULL) map.insert(std::pair<char,packet_type*>(type_code,type));
        else return map.find(type_code)->second;
    }
    template<class typeT>
    static void change_size(int field, int size){
        packet_type::get<typeT>()->change_field_size(field, size);
    }
    template<class typeT>
    static size_t get_size(int field = -1){
        if(field == -1) return packet_type::get<typeT>()->type_size;
        else return packet_type::get<typeT>()->sizes[field];
    }
    void change_field_size(int field, int size);
    void fill_packet(void* memory, const char type, va_list& fields);
    void commit();
protected:
    void construct(char code, int count, int* sizes, MPI_Datatype* types);
public:
    char type_code;
    int  type_size;                                // memory size of the type
    int  count;                                    // number of types used
    int* sizes;                                    // block-lengths of used types
    MPI_Aint* displacements;                       // memory placements of type-blocks
    MPI_Datatype mpi_type;                         // resulting MPI Datatype
    MPI_Datatype* compounds;                       // MPI compound datatypes
    BASE_FIELDS type;                              // mandatory packet-type first field
protected:
    packet_type();
};

#endif
