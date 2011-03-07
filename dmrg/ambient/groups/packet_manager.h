#ifndef AMBIENT_GROUPS_PACKET_MANAGER_H
#define AMBIENT_GROUPS_PACKET_MANAGER_H

#include "ambient/packets/types.h"
#include "ambient/packets/packet.h"

using namespace ambient::packets; 

namespace ambient{ namespace groups{

    class group;

    class ambient_request
    {
    public:
        ambient_request(void* memory);
        MPI_Request request;
        void* memory;
    };

    class packet_manager
    {
    private: 
        packet_manager(packet_manager const&);             // copy constructor is private
        packet_manager& operator=(packet_manager const&);  // assignment operator is private
        MPI_Comm* comm;
    public:
        std::vector<ambient_request*> send_requests;
        std::vector<ambient_request*> recv_requests;
        packet_manager(MPI_Comm* comm);
        void send(packet* pack, int dest);
        void recv(const packet_t& type, void* memory);
        ambient_request* isend(packet* pack, int dest);
        ambient_request* irecv(const packet_t& type, void* memory);
    };

    void wait_on_requests(group* grp);
    std::vector<ambient_request*>* get_send_requests(group* grp);
    std::vector<ambient_request*>* get_recv_requests(group* grp);

} }
#endif
