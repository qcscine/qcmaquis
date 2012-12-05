#ifndef AMBIENT_CHANNELS_MPI_PACKETS_DERIVED
#define AMBIENT_CHANNELS_MPI_PACKETS_DERIVED

// STANDARD AMBIENT FIELDS DEFINES
#define A_TYPE_FIELD 0 // MANDATORY FIRST TYPE CODE FIELD
#define A_DEST_FIELD 1 // RECOMMENDED SECOND FIELD IN DERIVED TYPES
#define A_OP_T_FIELD 2 // RECOMMENDED THIRD FIELD IN DERIVED TYPES
                       // (KEEP IT UNLESS YOU KNOW WHAT YOU ARE DOING)

// LAYOUT PACKET FIELDS DEFINES
#define A_LAYOUT_P_ACTION      3
#define A_LAYOUT_P_SID_FIELD   4
#define A_LAYOUT_P_STATE_FIELD 5
#define A_LAYOUT_P_OWNER_FIELD 6

#define A_BLOCK_P_SID_FIELD    3
#define A_BLOCK_P_STATE_FIELD  4
#define A_BLOCK_P_DATA_FIELD   5

#define A_CONTROL_P_SRC_FIELD  3
#define A_CONTROL_P_CODE_FIELD 4
#define A_CONTROL_P_INFO_FIELD 5

namespace ambient { namespace channels { namespace mpi {

// dest is mandatory first (in order to perform send operation
// w/o explicitely setting destination in send method). even if the 
// type is to be nested - dest is required to perform scatter op w/o 
// copying the data (if to do it implicitely).
    struct standard_packet_t : public packet_t
    {
        __a_fields__ dest, op_type;
        standard_packet_t(){
            __a_packet__
            dest     = MPI_INT;
            op_type  = MPI_BYTE;
            __a_pack{ 1, 1 }; 
            __a_code('1');
        }
    };

    struct control_packet_t : public standard_packet_t
    {
        __a_fields__ src, code, info;
        control_packet_t(){
            __a_packet__
            src      = MPI_INT;
            code     = MPI_BYTE;
            info     = MPI_BYTE;
            __a_pack{ 1, 1, 1 };
            __a_code('C');
        }
    };

    struct layout_packet_t : public standard_packet_t
    {
        __a_fields__ action, object_sid, state, owner;
        layout_packet_t(){
            __a_packet__
            action        = MPI_BYTE;
            object_sid    = MPI_INT;
            state         = MPI_BYTE;
            owner         = MPI_INT;
            __a_pack{ 1, 1, 1, 1 };
            __a_code('L');
        }
    };

    struct block_packet_t : public standard_packet_t
    {
        __a_flex_fields__ object_sid, state, data;
        block_packet_t(size_t size) 
        : standard_packet_t()
        {
            __a_packet__
            object_sid    = MPI_INT;
            state         = MPI_BYTE;
            data          = MPI_BYTE;
            __a_pack{ 1, 1, size };
            __a_code('Z'+size); // note: small size is dangerous (code overlap)
        }
    };

    class vbp { // volatile block packet (two-sided known communication)
    public:
        static vbp& instance(void* memory, size_t size, size_t dest){
            vbp& p = *(vbp*)memory;
            p.type = (double)((int)'Z' + (int)size/sizeof(double));
            p.dest = (double)dest;
            return p;
        }
        int length(){
            return (this->type-'Z');
        }
    public:
        double type;
        double dest;
        double data;
    };

} } }

#endif
