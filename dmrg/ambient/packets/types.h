#ifndef AMBIENT_PACKETS_TYPES_H
#define AMBIENT_PACKETS_TYPES_H

#include "ambient/packets/packet_t.h"

// STANDARD AMBIENT FIELDS DEFINES
#define A_TYPE_FIELD 0 // MANDATORY FIRST TYPE CODE FIELD
#define A_DEST_FIELD 1 // RECOMMENDED FIRST FIELD IN DERIVED TYPES
#define A_OP_T_FIELD 2 // RECOMMENDED SECOND FIELD IN DERIVED TYPES
                       // (KEEP IT UNLESS YOU KNOW WHAT YOU ARE DOING)

// LAYOUT PACKET FIELDS DEFINES
#define A_LAYOUT_P_OP_ID_FIELD 3
#define A_LAYOUT_P_OWNER_FIELD 4
#define A_LAYOUT_P_I_FIELD     5
#define A_LAYOUT_P_J_FIELD     6
#define A_LAYOUT_P_K_FIELD     7

#define A_BLOCK_P_OP_ID_FIELD  3
#define A_BLOCK_P_I_FIELD      4
#define A_BLOCK_P_J_FIELD      5
#define A_BLOCK_P_K_FIELD      6
#define A_BLOCK_P_DATA_FIELD   7

namespace ambient{ namespace packets{

// dest is mandatory first (in order to perform send operation
// w/o explicitely setting destination in send method). even if the 
// type is to be nested - dest is required to perform scatter op w/o 
// copying the data (if to do it implicitely).
    struct standard_packet_t : public packet_t
    {
        __a_fields__ dest, op_type;
        standard_packet_t()
        {
            __a_packet__
            dest     = MPI_INT;
            op_type  = MPI_BYTE;
            __a_pack{ 1, 1 }; 
            __a_code('1');
        }
    };

    struct control_packet_t : public standard_packet_t
    {
        __a_fields__ src, action, priority;
        control_packet_t()
        {
            __a_packet__
            src      = MPI_INT;
            action   = MPI_BYTE;
            priority = MPI_INT;
            __a_pack{ 1, 1, 1 };
            __a_code('C');
        }
    };

    struct layout_packet_t : public standard_packet_t
    {
        __a_fields__ op_profile_id, owner, i, j, k;
        layout_packet_t()
        {
            __a_packet__
            op_profile_id = MPI_INT;
            owner         = MPI_INT;
            i             = MPI_INT;
            j             = MPI_INT;
            k             = MPI_INT;
            __a_pack{ 1, 1, 1, 1, 1 };
            __a_code('L');
        }
    };

    struct block_packet_t : public standard_packet_t
    {
        __a_flex_fields__ op_profile_id, i, j, k, data;
        block_packet_t(size_t size) : standard_packet_t()
        {
            __a_packet__
            op_profile_id = MPI_INT;
            i             = MPI_INT;
            j             = MPI_INT;
            k             = MPI_INT;
            data          = MPI_DOUBLE;
            __a_pack{ 1, 1, 1, 1, size };
            __a_code('B');
        }
    };

} }
#endif
