#ifndef AMBIENT_PACKETS_H
#define AMBIENT_PACKETS_H

#include <memory.h>
#include <mpi.h>
#define LEN(x) sizeof(x)/sizeof(int)


class packet_type
{
private:
    packet_type(packet_type const&){};             // copy constructor is private
    packet_type& operator=(packet_type const&){};  // assignment operator is private
public:
    template<class typeT>
    static typeT* get(){
        static typeT* singleton = NULL;
        if(!singleton) singleton = new typeT();
        return singleton;
    }
    void instance(void* memory, ...);
    void change_field_size(int field, int size);

    int    size;                                   // memory size of the type
    int    count;                                  // number of types used
    int*   sizes;                                  // block-lengths of used types
    MPI_Aint* displacements;                       // memory placements of type-blocks
    MPI_Datatype mpi_type;                         // resulting MPI Datatype

    MPI_Datatype type, op_type, src, dst;          // packet fields
protected:
    packet_type(): type   (MPI_BYTE),              // this field should be always first!
                   op_type(MPI_BYTE),
                   src    (MPI_INT),
                   dst    (MPI_INT),
                   size(0), count(0), 
                   displacements(NULL), sizes(NULL)
    {
        static int sizes[] = { 1, 1, 1, 1 };
        packet_type::type_displacements(LEN(sizes), (MPI_Datatype*) &this->type, sizes);
    }
    void type_displacements(int count, MPI_Datatype* types, int* sizes);
};


class control_packet_type : public packet_type
{
    friend class packet_type;
    MPI_Datatype action, priority;
protected:
    control_packet_type(): action   (MPI_BYTE),
                           priority (MPI_INT)
    {
        static int sizes[] = { 1, 1 };
        packet_type::type_displacements(LEN(sizes), (MPI_Datatype*) &this->action, sizes);
    }
};

#endif
