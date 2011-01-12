#ifndef AMBIENT_PACKET_MANAGER_H
#define AMBIENT_PACKET_MANAGER_H

#include "ambient/packets/packet.h"
namespace ambient{ namespace packets{

    class packet_manager
    {
    private: 
        packet_manager();                                    // constructor is private
        packet_manager(packet_manager const&){};             // copy constructor is private
        packet_manager& operator=(packet_manager const&){};  // assignment operator is private
        MPI_Comm comm;
    public:
        static packet_manager* instance();
        void set_comm(MPI_Comm comm);
        void send(packet* pack, int dest);
        void recv(packet_t* type, void* memory);
    };

} }
#endif
