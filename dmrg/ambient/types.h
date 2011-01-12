#ifndef AMBIENT_TYPES_H
#define AMBIENT_TYPES_H

#include "ambient/packet_type.h"

// STANDARD AMBIENT FIELDS DEFINES
#define A_TYPE_FIELD 0 // MANDATORY FIRST TYPE CODE FIELD
#define A_DEST_FIELD 1 // RECOMMENDED FIRST FIELD IN DERIVED TYPES
                       // (KEEP IT UNLESS YOU KNOW WHAT YOU ARE DOING)


// dest is mandatory first (in order to perform send operation
// w/o explicitely setting destination in send method). even if the 
// type is to be nested - dest is required to perform scatter op w/o 
// copying the data (if to do it implicitely).
struct standard_packet_type : public packet_type
{
    FIELDS dest, op_type;
    standard_packet_type()
    {
        __A_PACKET__
        dest     = MPI_INT;
        op_type  = MPI_BYTE;
        PACK { 1, 1 }; 
        __A_CODE('1');
    }
};

struct control_packet_type : public standard_packet_type
{
    FIELDS src, action, priority;
    control_packet_type()
    {
        __A_PACKET__
        src      = MPI_INT;
        action   = MPI_BYTE;
        priority = MPI_INT;
        PACK {1, 1, 1};
        __A_CODE('C');
    }
};

struct data_packet_type : public standard_packet_type
{
    FIELDS src, priority, data;
    data_packet_type()
    {
        __A_PACKET__
        src      = MPI_INT;
        priority = MPI_INT;
        data     = MPI_DOUBLE;
        PACK {1, 1, 100};
        __A_CODE('D');
    }
};

#endif
