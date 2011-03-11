#include "ambient/ambient.h"
#include "ambient/packets/auxiliary.hpp"
#include "ambient/groups/packet_manager.h"
#include "ambient/groups/group.h"
#include "ambient/core/operation/operation.h"

namespace ambient{
#include "ambient/core/operation/operation.pp.sa.hpp"
}

namespace ambient{ namespace groups{

    void locking_handler(packet_manager::typed_q& in_q){
        ambient::packets::packet* pack = in_q.get_target_packet();
        char code = pack->get<char>(A_CONTROL_P_CODE_FIELD);
        if(code != 'L') return;
        char info = pack->get<char>(A_CONTROL_P_INFO_FIELD);
        if(in_q.manager->get_scope()->is_master()){
            if(info == 'C'){       // CLOSE
                in_q.manager->closure_mutex--;
                if(in_q.manager->closure_mutex == 0){
                    for(int i=0; i < in_q.manager->get_scope()->get_size(); i++)
                    in_q.manager->control_out->push(packets::pack<control_packet_t>(alloc_t<control_packet_t>(), i, "P2P",
                                                                                    ambient::rank(in_q.manager->get_scope()),
                                                                                    "LOCKING", "TRY TO CLOSE"));
                    in_q.manager->approve_closure_mutex = in_q.manager->get_scope()->get_size();
                }
            }else if(info == 'A'){ // APPROVE
                in_q.manager->approve_closure_mutex--;
                if(in_q.manager->approve_closure_mutex == 0){
                    for(int i=0; i < in_q.manager->get_scope()->get_size(); i++)
                    in_q.manager->control_out->push(packets::pack<control_packet_t>(alloc_t<control_packet_t>(), i, "P2P",
                                                                                    ambient::rank(in_q.manager->get_scope()),
                                                                                    "LOCKING", "FORCE CLOSURE"));
                    in_q.manager->closure_mutex = in_q.manager->get_scope()->get_size();
                }
            }else if(info == 'R'){ // REJECT
                    for(int i=0; i < in_q.manager->get_scope()->get_size(); i++)
                    in_q.manager->control_out->push(packets::pack<control_packet_t>(alloc_t<control_packet_t>(), i, "P2P",
                                                                                    ambient::rank(in_q.manager->get_scope()),
                                                                                    "LOCKING", "WITHDRAW CLOSURE"));
                    in_q.manager->closure_mutex = in_q.manager->get_scope()->get_size();
            }
        }
        if(info == 'T')      in_q.manager->state = packet_manager::CLOSURE; // TRY TO CLOSE
        else if(info == 'F') in_q.manager->state = packet_manager::CLOSED;  // FORCE CLOSURE
        else if(info == 'W') in_q.manager->state = packet_manager::OPEN;    // WITHDRAW CLOSURE
    }

    packet_manager::packet_manager(group* grp){
        this->scope = grp;
        this->comm = &grp->mpi_comm;
        this->control_in   = this->add_typed_q(get_t<control_packet_t>(), packet_manager::IN,  30);
        this->control_out  = this->add_typed_q(get_t<control_packet_t>(), packet_manager::OUT);
        this->layout_in_q  = this->add_typed_q(get_t<layout_packet_t>(),  packet_manager::IN, 30);
        this->layout_out_q = this->add_typed_q(get_t<layout_packet_t>(),  packet_manager::OUT);
        this->control_in->packet_delivered += new core::operation(locking_handler, this->control_in);
        this->state = packet_manager::OPEN;
        this->closure_mutex = this->scope->get_size();
    };
    void packet_manager::process(){
        size_t active_sends_number;
        for(;;){
            active_sends_number = 0;
            for(std::list<typed_q*>::const_iterator it = this->qs.begin(); it != qs.end(); ++it){
                if((*it)->active_requests_number > 0) for(int i=0; i < (*it)->priority; i++) (*it)->process();
                if((*it)->flow == OUT) active_sends_number += (*it)->active_requests_number;
            }

// finite state machine (closure proceedure).
// note: the state also can be modified in callback of control_in queue.
            if(this->state == OPEN && active_sends_number == 0){
                this->control_out->push(pack<control_packet_t>(alloc_t<control_packet_t>(), this->scope->get_master_g(), 
                                                               "P2P", ambient::rank(this->scope), 
                                                               "LOCKING", "CLOSURE")); 
                this->state = LOOSE;
            }else if(this->state == CLOSURE){ 
                if(active_sends_number == 0){
                    this->control_out->push(pack<control_packet_t>(alloc_t<control_packet_t>(), this->scope->get_master_g(), 
                                                                   "P2P", ambient::rank(this->scope),
                                                                   "LOCKING", "APPROVE")); 
                    this->state = LOOSE;
                }else{
                    this->control_out->push(pack<control_packet_t>(alloc_t<control_packet_t>(), this->scope->get_master_g(), 
                                                                   "P2P", ambient::rank(this->scope),
                                                                   "LOCKING", "REJECT")); 
                    this->state = OPEN;
                }
            }else if(this->state == CLOSED){
                this->state = OPEN; // leaving the door open ^_^
                printf("Closing the gate\n");
                return;
            }
        }
    }
    packet_manager::typed_q* packet_manager::add_typed_q(const packet_t& type, packet_manager::direction flow, int reservation, int priority){
        this->qs.push_back(new typed_q(this, type, flow, reservation, priority));
        return this->qs.back();
    }

    packet_manager::ambient_request::ambient_request(void* memory) 
    : memory(memory), request(MPI_REQUEST_NULL) { };

    packet_manager::typed_q::typed_q(packet_manager* manager, const packet_t& type, packet_manager::direction flow, int reservation, int priority) 
    : manager(manager), type(type), reservation(reservation), priority(priority), flow(flow), packet_delivered(), active_requests_number(0){
        if(flow == IN){ // make the initial recv request
            this->active_requests_number = this->reservation;
            for(int i=0; i < this->reservation; i++){
                this->requests.push_back(new ambient_request(alloc_t(this->type)));
                this->recv(this->requests.back());
            }
        }
    }

    void packet_manager::typed_q::push(packet* pack){
        this->active_requests_number++;
        assert(this->flow == OUT);
        for(int i=0; i < this->requests.size(); i++){ // locating free request
            if(this->requests[i]->request == MPI_REQUEST_NULL){
                delete((packet*)this->requests[i]->memory);
                this->requests[i]->memory = (void*)pack;
                this->send(this->requests[i]);
                return;
            }
        }
        this->requests.push_back(new ambient_request((void*)pack)); // subject for a change
        this->send(this->requests.back());
    }
    void packet_manager::typed_q::recv(ambient_request* request){
        MPI_Irecv(request->memory, 1, this->type.mpi_t, MPI_ANY_SOURCE, this->type.t_code, *this->manager->comm, &(request->request));
    }
    void packet_manager::typed_q::send(ambient_request* request){
        packet* pack = (packet*)request->memory;
        MPI_Isend(pack->data, 1, pack->mpi_t, *(int*)pack->get(A_DEST_FIELD), pack->get_t_code(), *this->manager->comm, &(request->request));
    }
    void packet_manager::typed_q::process(){
        int flag = 0;
        for(int i=0; i < this->requests.size(); i++)
        if(this->requests[i]->request != MPI_REQUEST_NULL){
            MPI_Test(&(this->requests[i]->request), &flag, MPI_STATUS_IGNORE);
            if(flag){
                if(this->flow == IN){ 
                    this->target_packet = unpack(this->type, this->requests[i]->memory); 
                    this->packet_delivered();
                    this->requests[i] = new ambient_request(alloc_t(this->type));
                    this->recv(this->requests[i]); // request renewal
//                    printf("Received info\n");
                }else if(this->flow == OUT){
                    this->active_requests_number--;
                    this->target_packet = (packet*)this->requests[i]->memory;
                    this->packet_delivered();
                }
            }
        }
    }

    packet* packet_manager::typed_q::get_target_packet(){
        return this->target_packet;
    }
    group* packet_manager::get_scope(){
        return this->scope;
    }
    packet_manager::typed_q::~typed_q(){ /* cancelling requests here */ }

} }
