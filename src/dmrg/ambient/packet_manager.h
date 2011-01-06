#ifndef AMBIENT_PACKET_MANAGER_H
#define AMBIENT_PACKET_MANAGER_H

#include "ambient/packets.h"

class packet_manager
{
private: 
    packet_manager();                                    // constructor is private
    packet_manager(packet_manager const&){};             // copy constructor is private
    packet_manager& operator=(packet_manager const&){};  // assignment operator is private
public:
    static packet_manager* instance();
    void commit(packet_type* type);
};


#endif
