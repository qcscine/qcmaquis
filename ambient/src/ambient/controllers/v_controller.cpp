#include "ambient/ambient.h"
#include "ambient/controllers/v_controller.h"
#include "ambient/channels/mpi_channel.h"
#include "ambient/channels/packets/auxiliary.hpp"
#include "ambient/utils/io.hpp"
#include "ambient/utils/timings.h"

#if __APPLE__ && __MACH__
#include <sched.h>	// for sched_yield()
#define pthread_yield() sched_yield()
#endif

#ifndef DEFAULT_NUM_THREADS 
#define DEFAULT_NUM_THREADS 1
#endif
#ifndef MAX_NUM_THREADS 
#define MAX_NUM_THREADS 8
#endif

// {{{ global objects accessible anywhere //
namespace ambient {
    channels::multirank& rank = channels::multirank::instance();
    models::v_model& model = models::v_model::instance();
    channels::ichannel& channel = channels::mpi_channel::instance();
    controllers::icontroller& controller = controllers::v_controller::instance();
    controllers::context& ctxt = controllers::context::instance();
    io cout;
    io cerr;
}
// }}} global objects accessible anywhere //

pthread_key_t pthread_env;
pthread_key_t pthread_tid;

namespace ambient { namespace controllers {

    v_controller::mod::mod(models::v_model::modifier* m, dim2 pin)
    : m(m), pin(pin) 
    {
    }

    v_controller::~v_controller(){
        for(size_t i = 1; i < this->num_threads; i++){
            this->tasks[i].active = false;
            pthread_join(this->pool[i], NULL);
            pthread_mutex_destroy(&this->mpool[i]);
        }
        pthread_mutex_destroy(&this->mutex);
        pthread_mutex_destroy(&this->pool_control_mutex);
    }

    v_controller::v_controller()
    : workload(0), rrn(0), num_threads(1) 
    {
        this->acquire(&ambient::channel);
        pthread_key_create(&pthread_env, free);
        pthread_key_create(&pthread_tid, free);
        pthread_mutex_init(&this->mutex, NULL);
        pthread_mutex_init(&this->pool_control_mutex, NULL);
        this->allocate_threads();
        this->set_num_threads(DEFAULT_NUM_THREADS);
    }

    pthread_mutex_t* v_controller::get_pool_control_mutex(){
        return &this->pool_control_mutex;
    }

    size_t v_controller::get_num_threads() const {
        return this->num_threads;
    }

    void v_controller::set_num_threads(size_t n){
        if(this->num_threads >= n || n > MAX_NUM_THREADS) return;
        for(size_t i = this->num_threads; i < n; i++){
            pthread_create(&this->pool[i], NULL, &v_controller::stream, &this->tasks[i]);
        }
        this->num_threads = n;
    }

    void v_controller::allocate_threads(){
        this->pool = (pthread_t*)malloc(sizeof(pthread_t)*MAX_NUM_THREADS);
        this->mpool = (pthread_mutex_t*)malloc(sizeof(pthread_mutex_t)*MAX_NUM_THREADS);
        this->tasks = new tasklist[MAX_NUM_THREADS];
        for(size_t i = 1; i < MAX_NUM_THREADS; i++){
            this->tasks[i].id = i;
            pthread_mutex_init(&this->mpool[i], NULL);
        }
        ctxt.set_tid(0); // master thread id is 0
    }

    void* v_controller::stream(void* list){
        mod* instruction;
        tasklist* l = static_cast<tasklist*>(list);
        ctxt.set_tid(l->id);

        while(l->active){
            instruction = (mod*)l->get_task();
            if(instruction == NULL){
                if(!l->idle){
                    l->idle = true;
                   /// printf("Thread is idle %lu\n", (size_t)l);
                }
                pthread_yield();
                continue;
            }else if(l->idle){
                l->idle = false;
                ///printf("Thread is not idle %lu\n", (size_t)l);
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

    void v_controller::push(models::v_model::modifier* op){
        this->stack.push_back(op);
        this->workload++; // should be done only in one thread
    }

    void v_controller::execute_mod(models::v_model::modifier* op, dim2 pin){
        if(op->pretend()) return;
        this->tasks[this->rrn].add_task(new mod(op, pin));
        ++this->rrn %= this->num_threads;
    }

    void v_controller::execute_free_mod(models::v_model::modifier* op){
        this->tasks[this->rrn].add_task(new mod(op, dim2(0,0)));
        ++this->rrn %= this->num_threads;
    }

    models::v_model::layout::entry* v_controller::alloc_block(models::v_model::revision& r){
        channels::packet_t& type = ambient::channel.get_block_packet_type(r.get_layout().get_mem_size());
        return new models::v_model::layout::entry(alloc_t(type), type.get_bound(A_BLOCK_P_DATA_FIELD));
    }

    models::v_model::layout::entry& v_controller::alloc_block(models::v_model::revision& r, size_t i, size_t j){
        channels::packet_t& type = ambient::channel.get_block_packet_type(r.get_layout().get_mem_size());
        r.get_layout().embed(alloc_t(type), i, j, type.get_bound(A_BLOCK_P_DATA_FIELD));
        return *r.block(i,j);
    }

    models::v_model::layout::entry& v_controller::ufetch_block(models::v_model::revision& r, size_t i, size_t j){
        if(r.block(i,j)->valid()){
            return *r.block(i,j);
        }else if(r.get_placement() == NULL || r.get_placement()->is_master()){
            assert(r.get_generator()->get_group() != NULL);
            if(r.get_generator()->get_group()->is_master()){
                this->alloc_block(r, i, j);
            }
        }else if(!r.block(i,j)->requested()){
            ambient::channel.ifetch(r.get_placement(), *r.get_layout().id().first, r.get_layout().id().second, i, j);
        }
        return *r.block(i,j);
    }

    models::v_model::layout::entry& v_controller::ifetch_block(models::v_model::revision& r, size_t i, size_t j){
        assert(r.get_placement() != NULL);
        if(r.block(i,j)->valid())
            this->atomic_receive(r, i, j); // check the stacked operations for the block
        else if(r.get_placement()->is_master()){
            if(r.get_layout().marked(i, j)){
                this->alloc_block(r, i, j);
                this->atomic_receive(r, i, j);
            }else{
                if(r.get_generator()->get_group() != NULL) // we already know the generation place
                    ambient::channel.ifetch(r.get_generator()->get_group(), *r.get_layout().id().first, r.get_layout().id().second, i, j);
//              else
//                  leaving on the side of the generator to deliver (we don't really know who generates but the destination is already known)
            }
        }else if(!r.block(i,j)->valid() && !r.block(i,j)->requested()){
            ambient::channel.ifetch(r.get_placement(), *r.get_layout().id().first, r.get_layout().id().second, i, j);
        }
        return *r.block(i,j);
    }

    void v_controller::atomic_receive(models::v_model::revision& r, size_t i, size_t j){
        // pthread_mutex_lock(&this->mutex); // will be needed in case of redunant accepts
        std::list<models::v_model::modifier*>::iterator it = r.block(i,j)->get_assignments().begin();
        while(it != r.block(i,j)->get_assignments().end()){
            this->execute_mod(*it, dim2(j,i));
            r.block(i,j)->get_assignments().erase(it++);
        }
    }

    void v_controller::atomic_complete(){
        pthread_mutex_lock(&this->mutex);
        assert(this->workload > 0);
        this->workload--;
        pthread_mutex_unlock(&this->mutex);
    }

    void v_controller::flush(){
        //static __a_timer time("ambient_total_compute_playout");
        //static __a_timer time2("ambient_total_logistic_playout");
        if(this->stack.empty()) return;
        while(!this->stack.end_reached())  // estimating operations credits 
            this->stack.pick()->weight();
        this->stack.sort();                // sorting operations using credit
        ctxt.state = context::EXECUTE;
        //time2.begin();
        while(!this->stack.end_reached()){
            ctxt.set_op(this->stack.pick());
            ctxt.get_op()->invoke();       // sending requests for data
        }
        //time2.end();
        //time.begin();
        this->master_stream(this->tasks);  // using up the main thread
        //time.end();
        ctxt.state = context::MARKUP;
        this->stack.clean();               // reseting the stack
    }
    
    channels::ichannel::packet* package(models::v_model::revision& r, const char* state, int i, int j, int dest){
        void* header = r.block(i,j)->get_memory();
        channels::ichannel::packet* package = channels::pack(channel.get_block_packet_type(r.get_layout().get_mem_size()), 
                                                             header, dest, "P2P", r.id().first, r.id().second, state, i, j, NULL);
        return package;
    }

    void forward_block(channels::ichannel::packet& cmd){
        channels::packet& c = static_cast<channels::packet&>(cmd);
        models::v_model::revision& r = *ambient::model.get_revision((size_t*)c.get(A_LAYOUT_P_GID_FIELD), 1, c.get<size_t>(A_LAYOUT_P_SID_FIELD));
        if(c.get<char>(A_LAYOUT_P_ACTION) != 'I') return; // INFORM OWNER ACTION
        size_t i = c.get<int>(A_LAYOUT_P_I_FIELD);
        size_t j = c.get<int>(A_LAYOUT_P_J_FIELD);
        models::v_model::layout::entry& entry = *r.block(i,j);
        if(entry.valid()){
            channel.emit(package(r, (const char*)c.get(A_LAYOUT_P_STATE_FIELD), i, j, c.get<int>(A_LAYOUT_P_OWNER_FIELD)));
        }else if(r.get_placement()->is_master() && r.get_layout().marked(i, j)){
            controller.alloc_block(r, i, j); // generating block
            forward_block(cmd);             // and forwarding
        }else{
            r.block(i,j)->get_path().push_back(c.get<int>(A_LAYOUT_P_OWNER_FIELD));
        }
    }

    void accept_block(channels::ichannel::packet& cmd){
        channels::packet& c = static_cast<channels::packet&>(cmd);
        size_t i = c.get<int>(A_BLOCK_P_I_FIELD);
        size_t j = c.get<int>(A_BLOCK_P_J_FIELD);
        models::v_model::revision& r = *ambient::model.get_revision((size_t*)c.get(A_BLOCK_P_GID_FIELD), 1, c.get<size_t>(A_BLOCK_P_SID_FIELD));
        if(r.block(i,j)->valid()) return; // quick exit for redunant accepts
        r.get_layout().embed(c.get_memory(), i, j, c.get_bound(A_BLOCK_P_DATA_FIELD));

        while(!r.block(i,j)->get_path().empty()){ // satisfying the path
            channel.emit(package(r, (const char*)c.get(A_LAYOUT_P_STATE_FIELD), 
                                 i, j, r.block(i,j)->get_path().back()));
            r.block(i,j)->get_path().pop_back();
        }

        controller.atomic_receive(r, i, j); // calling controller event handlers
    }

} }
