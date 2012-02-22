#include "ambient/ambient.h"
#include "ambient/controllers/v_controller.h"
#include "ambient/channels/mpi_channel.h"
#include "ambient/channels/packets/auxiliary.hpp"

#ifndef NUM_THREADS 
#define NUM_THREADS 4
#endif

// {{{ global objects accessible anywhere //
namespace ambient {
    channels::multirank& rank = channels::multirank::instance();
    models::imodel& model = models::v_model::instance();
    channels::ichannel& channel = channels::mpi_channel::instance();
    controllers::icontroller& controller = controllers::v_controller::instance();
    controllers::context& ctxt = controllers::context::instance();
}
// }}} global objects accessible anywhere //

pthread_key_t pthread_env;

namespace ambient { namespace controllers {

    v_controller::mod::mod(models::imodel::modifier* m, dim2 pin)
    : m(m), pin(pin) 
    {
    }

    v_controller::~v_controller(){
        for(int i = 1; i < NUM_THREADS; i++){
            this->tasks[i].active = false;
            pthread_join(this->pool[i], NULL);
        }
    }

    v_controller::v_controller()
    : workload(0) 
    {
        this->acquire(&ambient::channel);
        pthread_key_create(&pthread_env, free);
        pthread_mutex_init(&this->mutex, NULL);
        this->init_threads();
    }

    void v_controller::init_threads(){
        this->pool = (pthread_t*)malloc(sizeof(pthread_t)*NUM_THREADS);
        this->tasks = new tasklist[NUM_THREADS];
        for(int i = 1; i < NUM_THREADS; i++)
            pthread_create(&this->pool[i], NULL, &v_controller::stream, &this->tasks[i]);
    }

    void* v_controller::stream(void* list){
        mod* instruction;
        tasklist* l = static_cast<tasklist*>(list);

        while(l->active){
            instruction = (mod*)l->get_task();
            if(instruction == NULL){
                pthread_yield();
                continue;
            }
            ctxt.set_block_id(instruction->pin);
            instruction->m->invoke();
        }
        return NULL;
    }

    void v_controller::master_stream(void* list){
        mod* instruction;
        tasklist* l = static_cast<tasklist*>(list);

        while(this->workload){
            instruction = (mod*)l->get_task();
            if(instruction == NULL){
                pthread_yield();
                continue;
            }
            ctxt.set_block_id(instruction->pin);
            instruction->m->invoke();
        }
    }

    void v_controller::acquire(channels::ichannel* channel){
        this->channel = channel;
        this->channel->init();
    }

    void v_controller::push(models::imodel::modifier* op){
        this->stack.push_back(op);
        this->workload++;
    }

    void v_controller::push_mod(mod* m){
        static size_t rrn = 0; // ok since accessing only from MPI thread
        this->tasks[rrn].add_task(m);
        ++rrn %= NUM_THREADS;
    }

    models::imodel::layout::entry& v_controller::init_block(models::imodel::revision& r, size_t i, size_t j){
        channels::packet_t& type = ambient::channel.get_block_packet_type(r.get_layout().get_mem_size());
        r.get_layout().embed(alloc_t(type), i, j, type.get_bound(A_BLOCK_P_DATA_FIELD));
        
        ctxt.set_block_id(dim2(j,i)); // setting context block
        ((void(*)(models::v_model::object&))r.get_init())(static_cast<models::v_model::object&>(r.get_object()));
        return r(i,j);
    }

    models::imodel::layout::entry& v_controller::fetch_block(models::imodel::revision& r, size_t i, size_t j){
        assert(r.get_placement() != NULL);
        if(!r(i,j).valid() && !r(i,j).requested()) 
            ambient::channel.ifetch(r.get_placement(), *r.get_layout().id().first, r.get_layout().id().second, i, j);
        while(!r(i,j).valid()); // waiting for the receive to complete
        return r(i,j);
    }

    models::imodel::layout::entry& v_controller::ifetch_block(models::imodel::revision& r, size_t i, size_t j){
        if(r.get_placement() == NULL){
            this->init_block(r, i, j);
            this->atomic_receive(r, i, j);
        }else if(!r(i,j).valid() && !r(i,j).requested())
            ambient::channel.ifetch(r.get_placement(), *r.get_layout().id().first, r.get_layout().id().second, i, j);
        return r(i,j);
    }

    void v_controller::atomic_receive(models::imodel::revision& r, size_t i, size_t j){
        std::list<models::imodel::modifier*>::iterator it = r.get_modifiers().begin();
        while(it != r.get_modifiers().end()){
            if(&(*it)->get_pin() == &r)
                this->push_mod(new mod(*it, dim2(j,i)));
            it++;
            //r(i,j).get_assignments().erase(it++);
        }
    }

    void v_controller::atomic_complete(){
        pthread_mutex_lock(&this->mutex);
        this->workload--;
        pthread_mutex_unlock(&this->mutex);
    }

    void v_controller::flush(){
        if(this->stack.empty()) return;
        printf("Playout!\n");
        while(!this->stack.end_reached())  // estimating operations credits 
            this->stack.pick()->weight();
        this->stack.sort();                // sorting operations using credit
        ctxt.state = context::EXECUTE;
        printf("Executing logistics!\n");
        while(!this->stack.end_reached()){
            ctxt.set_op(this->stack.pick());
            ctxt.get_op()->invoke();       // sending requests for data
        }
        printf("Finished with the logistics!\n");
        this->master_stream(this->tasks);  // using up the main thread
        this->stack.clean();               // reseting the stack
    }
    
    channels::ichannel::packet* package(models::imodel::revision& r, const char* state, int i, int j, int dest){
        void* header = r(i,j).get_memory();
        channels::ichannel::packet* package = channels::pack(channel.get_block_packet_type(r.get_layout().get_mem_size()), 
                                                             header, dest, "P2P", r.id().first, r.id().second, state, i, j, NULL);
        return package;
    }

    void forward_block(channels::ichannel::packet& cmd){
        channels::packet& c = static_cast<channels::packet&>(cmd);
        models::imodel::revision& r = *ambient::model.get_revision((size_t*)c.get(A_LAYOUT_P_GID_FIELD), 1, c.get<size_t>(A_LAYOUT_P_SID_FIELD));
        if(c.get<char>(A_LAYOUT_P_ACTION) != 'I') return; // INFORM OWNER ACTION
        size_t i = c.get<int>(A_LAYOUT_P_I_FIELD);
        size_t j = c.get<int>(A_LAYOUT_P_J_FIELD);
        models::imodel::layout::entry& entry = r(i,j);
        if(entry.valid()){
            channel.emit(package(r, (const char*)c.get(A_LAYOUT_P_STATE_FIELD), i, j, c.get<int>(A_LAYOUT_P_OWNER_FIELD)));
        }else if(r.get_placement()->is_master() && r.get_layout().marked(i, j)){
            controller.init_block(r, i, j); // generating block
            forward_block(cmd);             // and forwarding
        }else{
            r.get_layout().push_path(i, j, c.get<int>(A_LAYOUT_P_OWNER_FIELD));
        }
    }

    void accept_block(channels::ichannel::packet& cmd){
        channels::packet& c = static_cast<channels::packet&>(cmd);
        size_t i = c.get<int>(A_BLOCK_P_I_FIELD);
        size_t j = c.get<int>(A_BLOCK_P_J_FIELD);
        models::imodel::revision& r = *ambient::model.get_revision((size_t*)c.get(A_BLOCK_P_GID_FIELD), 1, c.get<size_t>(A_BLOCK_P_SID_FIELD));
        if(r(i,j).valid()) return; // quick exit for redunant accepts
        r.get_layout().embed(c.get_memory(), i, j, c.get_bound(A_BLOCK_P_DATA_FIELD));

        while(!r.get_layout().get_path(i,j)->empty()){ // satisfying the path
            channel.emit(package(r, (const char*)c.get(A_LAYOUT_P_STATE_FIELD), 
                                 i, j, r.get_layout().pop_path(i,j)));
        }

        controller.atomic_receive(r, i, j); // calling controller event handlers
    }

} }
