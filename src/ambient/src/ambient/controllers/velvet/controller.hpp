#include "ambient/utils/io.hpp"
#include "ambient/utils/timings.hpp"

#if __APPLE__ && __MACH__
#include <sched.h>	// for sched_yield()
#define pthread_yield() sched_yield()
#endif

extern pthread_key_t pthread_tid;

namespace ambient { namespace controllers { namespace velvet {

    using ambient::channels::mpi::packet_t;
    using ambient::channels::mpi::packet;

    inline controller::~controller(){
        for(size_t i = 1; i < this->num_threads; i++){
            this->tasks[i].active = false;
            pthread_join(this->pool[i], NULL);
            pthread_mutex_destroy(&this->mpool[i]);
        }
        pthread_mutex_destroy(&this->mutex);
        pthread_mutex_destroy(&this->pool_control_mutex);
    }

    inline controller::controller()
    : workload(0), rrn(0), num_threads(1) 
    {
        this->acquire(&ambient::channel);
        pthread_key_create(&pthread_tid, free);
        pthread_mutex_init(&this->mutex, NULL);
        pthread_mutex_init(&this->pool_control_mutex, NULL);
        this->allocate_threads();
        this->set_num_threads(DEFAULT_NUM_THREADS);
    }

    inline pthread_mutex_t* controller::get_pool_control_mutex(){
        return &this->pool_control_mutex;
    }

    inline size_t controller::get_num_threads() const {
        return this->num_threads;
    }

    inline void controller::set_num_threads(size_t n){
        if(this->num_threads >= n || n > MAX_NUM_THREADS) return;
        for(size_t i = this->num_threads; i < n; i++){
            pthread_create(&this->pool[i], NULL, &controller::stream, &this->tasks[i]);
        }
        this->num_threads = n;
    }

    inline void controller::allocate_threads(){
        for(size_t i = 1; i < MAX_NUM_THREADS; i++){
            this->tasks[i].id = i;
            pthread_mutex_init(&this->mpool[i], NULL);
        }
        ctxt.set_tid(0); // master thread id is 0
    }

#ifndef AMBIENT_INTERFACE
    void* controller::stream(void* list){
        tasklist::task* instruction;
        tasklist* l = static_cast<tasklist*>(list);
        ctxt.set_tid(l->id);

        while(l->active){
            instruction = l->get_task();
            if(instruction == NULL){
                pthread_yield();
                continue;
            }
            ctxt.set_block_id(instruction->pin);
            instruction->f->computation();
            delete instruction;
        }
        return NULL;
    }
#endif

    inline void controller::master_stream(void* list){
        tasklist::task* instruction;
        tasklist* l = static_cast<tasklist*>(list);

        while(this->workload){
            instruction = l->get_task();
            if(instruction == NULL){
                printf("WARNING: MASTER HAS NULL INSTRUCTIONS!\n");
                pthread_yield();
                continue;
            }
            ctxt.set_block_id(instruction->pin);
            instruction->f->computation();
            delete instruction;
        }
    }

    inline void controller::acquire(channels::mpi::channel* channel){
        channel->init();
    }

    inline void controller::push(cfunctor* op){
        this->stack.push_back(op); // recasting to controllable
        this->workload++; // should be done only in one thread
    }

    inline void controller::execute_mod(cfunctor* op, dim2 pin){
        if(op->pretend()) return;
        this->tasks[this->rrn].add_task(new tasklist::task(op, pin));
        ++this->rrn %= this->num_threads;
    }

    inline void controller::execute_free_mod(cfunctor* op){
        this->tasks[this->rrn].add_task(new tasklist::task(op, dim2(0,0)));
        ++this->rrn %= this->num_threads;
    }

    inline void controller::alloc_block(const memspec& s, layout& l, size_t x, size_t y){
        l.embed(s.alloc(), x, y, s.get_bound());
    }

    inline void controller::calloc_block(const memspec& s, layout& l, size_t x, size_t y){
        l.embed(s.calloc(), x, y, s.get_bound());
    }

    // note: ufetch_block is used only by pt_fetch in user-space
    inline layout::entry& controller::ufetch_block(revision& r, size_t x, size_t y){
        //if(r.block(x,y).valid()){
            return r.block(x,y);
        //}else printf("REQUESTING UNEXISTING BLOCK (%lu, %lu)!\n", x, y);
            
            //else if(r.get_placement() == NULL || r.get_placement()->is_master()){
            //assert(r.get_generator()->get_group() != NULL);
            //if(r.get_generator()->get_group()->is_master()){
            //  return this->alloc_block(r.get_layout(), x, y);
            //}
        //}
        //else if(!r.block(x,y).requested()){
         //   ambient::channel.ifetch(r.get_placement(), *r.get_layout().id(), x, y);
         //} // blocks should be already requested
        //return r.block(x,y);
    }

    inline layout::entry& controller::ifetch_block(revision& r, size_t x, size_t y){
        this->atomic_receive(r.get_layout(), x, y); // check the stacked operations for the block
        return r.block(x,y);

        /*assert(r.get_placement() != NULL);
        if(r.block(x,y).valid())
            this->atomic_receive(r.get_layout(), x, y); // check the stacked operations for the block
        else if(r.get_placement()->is_master()){
            if(r.get_layout().marked(x, y)){
                this->alloc_block(r.get_layout(), x, y); // ! making all allocations inside computational kernels
                this->atomic_receive(r.get_layout(), x, y);
            }else{
                printf("UNEXPECTED ERROR IN IFETCH BLOCK -- asking for %lu %lu!\n", x, y);
                if(r.get_generator()->get_group() != NULL) // we already know the generation place
                    ambient::channel.ifetch(r.get_generator()->get_group(), r.get_layout().id(), x, y);
              //else
              //    leaving on the side of the generator to deliver (we don't really know who generates but the destination is already known)
            }
        }else if(!r.block(x,y).valid() && !r.block(x,y).requested()){
            ambient::channel.ifetch(r.get_placement(), r.get_layout().id(), x, y);
        }
        return r.block(x,y);*/
    }

    inline bool controller::lock_block(revision& r, size_t x, size_t y){
#ifdef AMBIENT_THREADS
        pthread_mutex_lock(&this->pool_control_mutex);
        bool acquired = r.block(x,y).trylock();
        pthread_mutex_unlock(&this->pool_control_mutex);
        return acquired;
#else
        return true;
#endif
    }

    inline void controller::unlock_block(revision& r, size_t x, size_t y){
#ifdef AMBIENT_THREADS
        pthread_mutex_lock(&this->pool_control_mutex);
        r.block(x,y).unlock();
        pthread_mutex_unlock(&this->pool_control_mutex);
#endif
    }

    inline void controller::atomic_receive(layout& l, size_t x, size_t y){
        // pthread_mutex_lock(&this->mutex); // will be needed in case of redunant accepts
        std::list<cfunctor*>::iterator it = l.get(x,y).get_assignments().begin();
        while(it != l.get(x,y).get_assignments().end()){
            this->execute_mod(*it, dim2(x,y));
            l.get(x,y).get_assignments().erase(it++);
        }
    }

    inline void controller::atomic_complete(cfunctor* op){
#ifdef AMBIENT_THREADS
        pthread_mutex_lock(&this->mutex);
#endif
        assert(this->workload > 0);
        this->workload--;
#ifdef AMBIENT_THREADS
        pthread_mutex_unlock(&this->mutex);
#endif
        delete op;
    }

    inline void controller::conditional_flush(){
        //static int counter = 1;
        //counter = (counter+1) % 2;
        //if(!counter) 
            this->flush();
    }

    inline void controller::flush(){
        //static __a_timer time("ambient_total_compute_playout");
        //static __a_timer time2("ambient_total_logistic_playout");
        //if(this->stack.empty()) return;

        //while(!this->stack.end_reached())  // estimating operations credits 
        //    this->stack.pick()->weight();
        //this->stack.sort();                // sorting operations using credit
        //time2.begin();
        while(!this->stack.end_reached())    // can be read only once
            this->stack.pick()->logistics(); // sending requests for data
        //time2.end();
        //time.begin();
        this->master_stream(this->tasks);  // using up the main thread
        this->stack.reset();
        //time.end();
    }
    
    inline packet* package(layout& l, const char* state, int x, int y, int dest){
        void* header = l.get(x,y).get_memory();
        //if(header == NULL) printf("HEADER IS NULL (SWAPPED)\n");
        packet* package = pack(*(packet_t*)l.get_spec().get_packet_t(), 
                               header, dest, "P2P", l.sid, state, x, y, NULL);
        return package;
    }

    inline void forward_block(packet& cmd){
        packet& c = static_cast<packet&>(cmd);
        layout& l = *ambient::model.get_layout(c.get<size_t>(A_LAYOUT_P_SID_FIELD));
        if(c.get<char>(A_LAYOUT_P_ACTION) != 'I') return; // INFORM OWNER ACTION
        size_t x = c.get<int>(A_LAYOUT_P_X_FIELD);
        size_t y = c.get<int>(A_LAYOUT_P_Y_FIELD);
        layout::entry& entry = l.get(x,y);
        if(entry.valid()){
            channel.emit(package(l, (const char*)c.get(A_LAYOUT_P_STATE_FIELD), x, y, c.get<int>(A_LAYOUT_P_OWNER_FIELD)));
        }else if(l.placement->is_master()){
            ambient::controller.alloc_block(l.spec, l, x, y); // generating block
            forward_block(cmd);             // and forwarding
        }else{
            l.get(x,y).get_path().push_back(c.get<int>(A_LAYOUT_P_OWNER_FIELD));
        }
    }

    inline void accept_block(packet& cmd){
        packet& c = static_cast<packet&>(cmd);
        size_t x = c.get<int>(A_BLOCK_P_X_FIELD);
        size_t y = c.get<int>(A_BLOCK_P_Y_FIELD);
        layout& l = *ambient::model.get_layout(c.get<size_t>(A_BLOCK_P_SID_FIELD));
        if(l.get(x,y).valid()) return; // quick exit for redunant accepts
        l.embed(c.get_memory(), x, y, c.get_bound(A_BLOCK_P_DATA_FIELD));

        while(!l.get(x,y).get_path().empty()){ // satisfying the path
            channel.emit(package(l, (const char*)c.get(A_LAYOUT_P_STATE_FIELD), 
                                 x, y, l.get(x,y).get_path().back()));
            l.get(x,y).get_path().pop_back();
        }

        ambient::controller.atomic_receive(l, x, y); // calling controller event handlers
    }

} } }
