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
    public:
        static packet_manager* instance();
        void send(packet* pack, int dest);
    };

} }
#endif
