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
                this->emit(pack<control_packet_t>(alloc_t<control_packet_t>(), this->grp->get_master(), 
                                                  "P2P", ambient::rank(this->grp), "LOCKING", "CLOSURE")); 
                this->state = LOOSE;
        }else if(this->state == CLOSURE){ 
            if(active_sends_number == 0){
                this->emit(pack<control_packet_t>(alloc_t<control_packet_t>(), this->grp->get_master(), 
                                                  "P2P", ambient::rank(this->grp), "LOCKING", "APPROVE")); 
                this->state = LOOSE;
            }else{
                this->emit(pack<control_packet_t>(alloc_t<control_packet_t>(), this->grp->get_master(), 
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
        this->add_typed_q(type, packet_manager::IN,  60);
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
            ambient::spin();
            active_sends_number = 0;
            this->spin();
            for(std::list<typed_q*>::const_iterator it = this->qs.begin(); it != qs.end(); ++it)
                if((*it)->flow == OUT) active_sends_number += (*it)->get_active();
            if(this->process_locking(active_sends_number)){ return this->spin(); } // spinning last time
            //counter++;
            //if(counter == 1000) printf("R%d: inactive for %d iterations... (active sends: %d)\n", ambient::rank(), (int)counter, (int)active_sends_number);
        }
    }
    void packet_manager::spin(int n){
        for(int i=0; i < n; i++){
            for(std::list<typed_q*>::const_iterator it = this->qs.begin(); it != qs.end(); ++it){
                if((*it)->get_active() > 0) for(int i=0; i < (*it)->priority; i++) (*it)->spin();
            }
        }
    }
    packet_manager::typed_q* packet_manager::add_typed_q(const packet_t& type, packet_manager::direction flow, int reservation, int priority){
        this->qs.push_back(new typed_q(this, type, flow, reservation, priority));
        return this->qs.back();
    }

    packet_manager::request::request(void* memory) 
    : memory(memory), mpi_request(MPI_REQUEST_NULL), fail_count(0) { };

    packet_manager::request* packet_manager::typed_q::get_request(){
        if(this->free_reqs.size() == 0){
            this->free_reqs.push_back(new request(alloc_t(this->type)));
        }
        request* handle = this->free_reqs.back();
        this->free_reqs.pop_back();
        return handle;
    }
    void packet_manager::typed_q::return_request(request* r){
        this->free_reqs.push_back(r);
    }
    packet_manager::typed_q::typed_q(packet_manager* manager, const packet_t& type, packet_manager::direction flow, int reservation, int priority) 
    : manager(manager), type(type), priority(priority), flow(flow), packet_delivered()
    {    
        for(int i=0; i < PULL_RESERVATION; i++){
            this->free_reqs.push_back(new request(alloc_t(this->type)));
        }
        if(flow == IN){ // make the initial recv request
            for(int i=0; i < reservation; i++){
                this->reqs.push_back(this->get_request());
                this->recv(this->reqs.back());
            }
        }
    }
    void packet_manager::typed_q::push(packet* pack){
        assert(this->flow == OUT);
        this->reqs.push_back(this->get_request());
        request* target = this->reqs.back();
        target->memory = (void*)pack; // don't forget to free memory (if not needed)
        this->send(target);
    }
    void packet_manager::typed_q::recv(request* r){
        MPI_Irecv(r->memory, 1, this->type.mpi_t, MPI_ANY_SOURCE, this->type.t_code, *this->manager->comm, &(r->mpi_request));
    }
    void packet_manager::typed_q::send(request* r){
        packet* pack = (packet*)r->memory;
        if(pack->get<char>(A_OP_T_FIELD) == 'P'){
            assert(pack->get<int>(A_DEST_FIELD) >= 0);
            MPI_Isend(pack->data, 1, pack->mpi_t, pack->get<int>(A_DEST_FIELD), pack->get_t_code(), *this->manager->comm, &(r->mpi_request));
        }else if(pack->get<char>(A_OP_T_FIELD) == 'B'){
            pack->set(A_OP_T_FIELD, "P2P");
            pack->set(A_DEST_FIELD, 0); // using up this request
            MPI_Isend(pack->data, 1, pack->mpi_t, pack->get<int>(A_DEST_FIELD), pack->get_t_code(), *this->manager->comm, &(r->mpi_request));
            for(size_t i=1; i < this->manager->get_group()->get_size(); i++){
                pack->set(A_DEST_FIELD, i);
                this->push(pack);
            }
        }
    }
    size_t packet_manager::typed_q::get_active(){
        if(this->flow == IN) return 1;
        return this->reqs.size();
    }
    void packet_manager::typed_q::spin(){
        int flag = 0;
        std::list<request*>::iterator it;
        if(this->flow == IN)
        for(it = this->reqs.begin(); it != this->reqs.end(); ++it){
            MPI_Test(&((*it)->mpi_request), &flag, MPI_STATUS_IGNORE);
            if(flag){
                this->target_packet = unpack(this->type, (*it)->memory); 
                this->packet_delivered();
                (*it)->memory = alloc_t(this->type); // otherwise will overwrite old memory
                this->recv(*it); // request renewal
            }
        }
        else 
        for(it = this->reqs.begin(); it != this->reqs.end(); ++it){
            MPI_Test(&((*it)->mpi_request), &flag, MPI_STATUS_IGNORE);
            if(flag){
                this->target_packet = (packet*)(*it)->memory;
                this->packet_delivered();
                this->return_request(*it);
                it = this->reqs.erase(it);
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
