#include "ambient/ambient.h"
#include "ambient/packets/auxiliary.hpp"
#include "ambient/groups/packet_manager.h"
#include "ambient/groups/group.h"
#include "ambient/core/operation/operation.h"
#include "ambient/core/operation/operation.pp.sa.hpp"

namespace ambient{ namespace core{
using namespace ambient::groups;

    extern void update_layout(packet_manager::typed_q& in_q);
    extern void forward_layout(packet_manager::typed_q& in_q);
    extern void forward_block(packet_manager::typed_q& in_q);
} }

namespace ambient{ namespace groups{

    void locking_handler(packet_manager::typed_q& in_q){
        ambient::packets::packet* pack = in_q.get_target_packet();
        char code = pack->get<char>(A_CONTROL_P_CODE_FIELD);
        if(code != 'L') return;
        char info = pack->get<char>(A_CONTROL_P_INFO_FIELD);
        if(in_q.manager->get_group()->is_master()){
            if(info == 'C'){       // CLOSE
                in_q.manager->closure_mutex--;
                if(in_q.manager->closure_mutex == 0){
                    for(int i=0; i < in_q.manager->get_group()->get_size(); i++)
                    in_q.manager->emit(packets::pack<control_packet_t>(alloc_t<control_packet_t>(), i, "P2P",
                                                                       ambient::rank(in_q.manager->get_group()),
                                                                       "LOCKING", "TRY TO CLOSE"));
                    in_q.manager->approve_closure_mutex = in_q.manager->get_group()->get_size();
                }
            }else if(info == 'A'){ // APPROVE
                in_q.manager->approve_closure_mutex--;
                if(in_q.manager->approve_closure_mutex == 0){
                    for(int i=0; i < in_q.manager->get_group()->get_size(); i++)
                    in_q.manager->emit(packets::pack<control_packet_t>(alloc_t<control_packet_t>(), i, "P2P",
                                                                       ambient::rank(in_q.manager->get_group()),
                                                                       "LOCKING", "FORCE CLOSURE"));
                    in_q.manager->closure_mutex = in_q.manager->get_group()->get_size();
                }
            }else if(info == 'R'){ // REJECT
                    for(int i=0; i < in_q.manager->get_group()->get_size(); i++)
                    in_q.manager->emit(packets::pack<control_packet_t>(alloc_t<control_packet_t>(), i, "P2P",
                                                                       ambient::rank(in_q.manager->get_group()),
                                                                       "LOCKING", "WITHDRAW CLOSURE"));
                    in_q.manager->closure_mutex = in_q.manager->get_group()->get_size();
            }
        }
        if(info == 'T')      in_q.manager->state = packet_manager::CLOSURE; // TRY TO CLOSE
        else if(info == 'F') in_q.manager->state = packet_manager::CLOSED;  // FORCE CLOSURE
        else if(info == 'W') in_q.manager->state = packet_manager::OPEN;    // WITHDRAW CLOSURE
    }

    bool packet_manager::process_locking(size_t active_sends_number){
// finite state machine (closure proceedure).
// note: the state also can be modified in callback of control_in queue.
        if(this->state == OPEN && active_sends_number == 0){
                this->emit(pack<control_packet_t>(alloc_t<control_packet_t>(), this->grp->get_master_g(), 
                                                  "P2P", ambient::rank(this->grp), "LOCKING", "CLOSURE")); 
                this->state = LOOSE;
        }else if(this->state == CLOSURE){ 
            if(active_sends_number == 0){
                this->emit(pack<control_packet_t>(alloc_t<control_packet_t>(), this->grp->get_master_g(), 
                                                  "P2P", ambient::rank(this->grp), "LOCKING", "APPROVE")); 
                this->state = LOOSE;
            }else{
                this->emit(pack<control_packet_t>(alloc_t<control_packet_t>(), this->grp->get_master_g(), 
                                                  "P2P", ambient::rank(this->grp), "LOCKING", "REJECT")); 
                this->state = OPEN;
            }
        }else if(this->state == CLOSED){
            this->state = OPEN; // leaving the door open ^_^
            return true;
        }
        return false;
    }

    bool packet_manager::subscribed(const packet_t& type){
        for(std::list<typed_q*>::const_iterator it = this->qs.begin(); it != this->qs.end(); ++it){
            if(&((*it)->type) == &type) return true;
        }
        return false;
    }

    void packet_manager::subscribe(const packet_t& type){
        this->add_typed_q(type, packet_manager::IN,  30);
    }

    void packet_manager::add_handler(const packet_t& type, core::operation* callback){
        for(std::list<typed_q*>::const_iterator it = this->qs.begin(); it != qs.end(); ++it){
            if(&((*it)->type) == &type && (*it)->flow == IN) (*it)->packet_delivered += callback;
        }
    }

    void packet_manager::emit(packet* pack){
        for(std::list<typed_q*>::const_iterator it = this->qs.begin(); it != qs.end(); ++it){
            if(&((*it)->type) == &(pack->get_t()) && (*it)->flow == OUT){ (*it)->push(pack); return; }
        }
        this->add_typed_q(pack->get_t(), packet_manager::OUT)->push(pack);
    }

    packet_manager::typed_q* packet_manager::get_pipe(const packet_t& type, packet_manager::direction flow){
        for(std::list<typed_q*>::const_iterator it = this->qs.begin(); it != qs.end(); ++it){
            if(&((*it)->type) == &(type) && flow == (*it)->flow){ return (*it); }
        }
        return this->add_typed_q(type, flow);
    }

    packet_manager::packet_manager(group* grp) {
        this->grp = grp;
        this->comm = &grp->mpi_comm;

        this->subscribe(get_t<control_packet_t>());
        this->subscribe(get_t<layout_packet_t>() );
        this->add_handler( get_t<control_packet_t>(), new core::operation(locking_handler     , this->get_pipe(get_t<control_packet_t>(), IN)) );
        this->add_handler( get_t<layout_packet_t>() , new core::operation(core::forward_block , this->get_pipe(get_t<layout_packet_t>() , IN)) );
        this->add_handler( get_t<layout_packet_t>() , new core::operation(core::forward_layout, this->get_pipe(get_t<layout_packet_t>() , IN)) );
        this->add_handler( get_t<layout_packet_t>() , new core::operation(core::update_layout , this->get_pipe(get_t<layout_packet_t>() , IN)) );
// note: the order of handlers really matters
        this->state = packet_manager::OPEN;
        this->closure_mutex = this->grp->get_size();
    };
    void packet_manager::spin_loop(){
        size_t active_sends_number;
        size_t counter = 0;
        for(;;){
            active_sends_number = 0;
            for(std::list<typed_q*>::const_iterator it = this->qs.begin(); it != qs.end(); ++it){
                if((*it)->active_requests_number > 0) for(int i=0; i < (*it)->priority; i++) (*it)->spin();
                if((*it)->flow == OUT) active_sends_number += (*it)->active_requests_number;
            }
            if(this->process_locking(active_sends_number)){ this->spin(); break; } // spinning last time
            counter++;
            if(counter == 1000) printf("R%d: I'm stuck in first!!! (active sends: %d)\n", ambient::rank(), active_sends_number);
        }
        for(;;){
            active_sends_number = 0;
            for(std::list<typed_q*>::const_iterator it = this->qs.begin(); it != qs.end(); ++it){
                if((*it)->active_requests_number > 0) for(int i=0; i < (*it)->priority; i++) (*it)->spin();
                if((*it)->flow == OUT) active_sends_number += (*it)->active_requests_number;
            }
            if(this->process_locking(active_sends_number)){ this->spin(); break; } // spinning last time
            counter++;
            if(counter == 1000) printf("R%d: I'm stuck in second!!! (active sends: %d)\n", ambient::rank(), active_sends_number);
        }
        MPI_Barrier(*this->comm);
    }
    void packet_manager::spin(int n){
        for(int i=0; i < n; i++){
            for(std::list<typed_q*>::const_iterator it = this->qs.begin(); it != qs.end(); ++it){
                if((*it)->active_requests_number > 0) for(int i=0; i < (*it)->priority; i++) (*it)->spin();
            }
        }
    }
    packet_manager::typed_q* packet_manager::add_typed_q(const packet_t& type, packet_manager::direction flow, int reservation, int priority){
        this->qs.push_back(new typed_q(this, type, flow, reservation, priority));
        return this->qs.back();
    }

    packet_manager::ambient_request::ambient_request(void* memory) 
    : memory(memory), request(MPI_REQUEST_NULL), fail_count(0) { };

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
                this->requests[i]->memory = (void*)pack; // don't forget to free memory (if not needed)
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
        if(pack->get<char>(A_OP_T_FIELD) == 'P'){
            MPI_Isend(pack->data, 1, pack->mpi_t, *(int*)pack->get(A_DEST_FIELD), pack->get_t_code(), *this->manager->comm, &(request->request));
        }else if(pack->get<char>(A_OP_T_FIELD) == 'B'){
            pack->set(A_OP_T_FIELD, "P2P");
            pack->set(A_DEST_FIELD, 0); // using up this request
            MPI_Isend(pack->data, 1, pack->mpi_t, *(int*)pack->get(A_DEST_FIELD), pack->get_t_code(), *this->manager->comm, &(request->request));
            for(int i=1; i < this->manager->get_group()->get_size(); i++){
                pack->set(A_DEST_FIELD, i);
                this->push(pack);
            }
        }
    }
    void packet_manager::typed_q::spin(){
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
                }else if(this->flow == OUT){
                    this->active_requests_number--;
                    this->target_packet = (packet*)this->requests[i]->memory;
                    this->packet_delivered();
                    this->requests[i]->fail_count = 0;
                }
            }else{
                if(this->flow == OUT){ 
                    this->requests[i]->fail_count++;
                    if(this->requests[i]->fail_count == 200) printf("The failed request: %d\n", i);
                }
            }
        }
    }

    packet* packet_manager::typed_q::get_target_packet(){
        return this->target_packet;
    }
    group* packet_manager::get_group(){
        return this->grp;
    }
    packet_manager::typed_q::~typed_q(){ /* cancelling requests here */ }

} }
