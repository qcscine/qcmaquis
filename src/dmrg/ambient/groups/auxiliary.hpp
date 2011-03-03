#ifndef AMBIENT_GROUPS_AUX_H
#define AMBIENT_GROUPS_AUX_H

namespace ambient{ namespace groups{

    template<typename T>
    packet* recv(group* grp, void* memory){
        grp->manager->recv(get_t<T>(), memory);
        return unpack<T>(memory);
    }

    template<typename T>
    packet* recv(const char* grp, void* memory)
    {
        return recv<T>(group_map(grp), memory);
    }

    packet* recv(const packet_t& type, group* grp, void* memory);
    packet* recv(const packet_t& type, const char* grp, void* memory);
    ambient_request* irecv(const packet_t& type, group* grp, void* memory);
    ambient_request* irecv(const packet_t& type, const char* grp, void* memory);
    void send(packet* pack, group* grp, int dest = -1);
    void send(packet* pack, const char* grp, int dest = -1);
    ambient_request* isend(packet* pack, group* grp, int dest = -1);
    ambient_request* isend(packet* pack, const char* grp, int dest = -1);

} }

#endif
