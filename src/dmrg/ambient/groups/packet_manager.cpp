#include "ambient/ambient.h"
#include "ambient/groups/packet_manager.h"
#include "ambient/groups/group.h"

namespace ambient{ namespace groups{

    ambient_request::ambient_request(void* memory) : memory(memory), request(MPI_REQUEST_NULL) {};

    packet_manager::packet_manager(MPI_Comm* comm):send_requests(),recv_requests(){
        this->comm = comm;
    };

    void packet_manager::send(packet* pack, int dest)
    {
        if(pack->get_mpi_t() != pack->mpi_t) printf("Erroneous packet type: the type has changed since packet's creation\n");
        MPI_Send(pack->data, 1, pack->mpi_t, dest, pack->get_t_code(), *this->comm);
    }

    ambient_request* packet_manager::isend(packet* pack, int dest)
    {
        if(pack->get_mpi_t() != pack->mpi_t) printf("Erroneous packet type: the type has changed since packet's creation\n");
        this->send_requests.push_back(new ambient_request((void*)pack));
        MPI_Isend(pack->data, 1, pack->mpi_t, dest, pack->get_t_code(), *this->comm, &(send_requests.back()->request));
        return send_requests.back();
    }

    void packet_manager::recv(const packet_t& type, void* memory)
    {
        MPI_Recv(memory, 1, type.mpi_t, MPI_ANY_SOURCE, type.t_code, *this->comm, MPI_STATUS_IGNORE);
    }

    ambient_request* packet_manager::irecv(const packet_t& type, void* memory)
    {
        this->recv_requests.push_back(new ambient_request(memory));
        MPI_Irecv(memory, 1, type.mpi_t, MPI_ANY_SOURCE, type.t_code, *this->comm, &(recv_requests.back()->request));
        return recv_requests.back();
    }

    void send(packet* pack, group* grp, int dest = -1)
    {
        if(dest != -1){
            grp->manager->send(pack, dest);
        }else{
            if(pack->get_t().compounds[1] != MPI_INT) printf("Warning: the dest field (#1) is not of type MPI_INT!\n");
            grp->manager->send(pack, *(int*)pack->get(A_DEST_FIELD));
        }
    }
    ambient_request* isend(packet* pack, group* grp, int dest = -1)
    {
        if(dest != -1) return grp->manager->isend(pack, dest);
        if(pack->get_t().compounds[1] != MPI_INT) printf("Warning: the dest field (#1) is not of type MPI_INT!\n");
        return grp->manager->isend(pack, *(int*)pack->get(A_DEST_FIELD));
    }

    packet* recv(const packet_t& type, group* grp, void* memory){
        grp->manager->recv(type, memory);
        return unpack(type, memory);
    }
    ambient_request* irecv(const packet_t& type, group* grp, void* memory){
        return grp->manager->irecv(type, memory);
    }

    void wait_on_requests(group* grp){
        int flag = 0;
        int ir = 0;
        int is = 0;
        int counter = 0;
        int count_recv = grp->manager->recv_requests.size();
        int count_send = grp->manager->send_requests.size();
        int goal = count_recv + count_send;

        while(counter != goal){
            if(count_send){
                if(grp->manager->send_requests[is]->request != MPI_REQUEST_NULL){
                    MPI_Test(&(grp->manager->send_requests[is]->request), &flag, MPI_STATUS_IGNORE);
                    if(flag) counter++;
                }
                is = (is+1) % count_send;
            }
            if(count_recv){
                if(grp->manager->recv_requests[ir]->request != MPI_REQUEST_NULL){
                    MPI_Test(&(grp->manager->recv_requests[ir]->request), &flag, MPI_STATUS_IGNORE);
                    if(flag) counter++;
                }
                ir = (ir+1) % count_recv;
            }
        }
    }

    std::vector<ambient_request*>* get_send_requests(group* grp){
        return &grp->manager->send_requests;
    }

    std::vector<ambient_request*>* get_recv_requests(group* grp){
        return &grp->manager->recv_requests;
    }

    void send(packet* pack, const char* grp, int dest = -1){  send(pack, group_map(grp), dest); }
    packet* recv(const packet_t& type, const char* grp, void* memory){ return recv(type, group_map(grp), memory); }
    ambient_request* isend(packet* pack, const char* grp, int dest = -1){ isend(pack, group_map(grp), dest); }
    ambient_request* irecv(const packet_t& type, const char* grp, void* memory){ return irecv(type, group_map(grp), memory); }
} }
