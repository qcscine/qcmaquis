#ifndef AMBIENT_GROUPS_PACKET_MANAGER_H
#define AMBIENT_GROUPS_PACKET_MANAGER_H

#include "ambient/packets/types.h"
#include "ambient/packets/packet.h"

using namespace ambient::packets; 

namespace ambient{ namespace groups{

    class packet_manager
    {
    private: 
        packet_manager(packet_manager const&){};             // copy constructor is private
        packet_manager& operator=(packet_manager const&){};  // assignment operator is private
        MPI_Comm* comm;
    public:
        packet_manager(MPI_Comm* comm);
        void send(packet* pack, int dest);
        void recv(const packet_t& type, void* memory);
    };

} }
#endif
