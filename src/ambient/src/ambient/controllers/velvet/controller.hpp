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
        }
        pthread_mutex_destroy(&this->mutex);
    }

    inline controller::controller()
    : workload(0), rrn(0), num_threads(1), muted(false) 
    {
        this->acquire(&ambient::channel);
        pthread_key_create(&pthread_tid, free);
        pthread_mutex_init(&this->mutex, NULL);
        this->allocate_threads();
        this->set_num_threads(AMBIENT_THREADS);
    }

    inline size_t controller::get_num_threads() const {
        return this->num_threads;
    }

    inline void controller::set_num_threads(size_t n){
        if(this->num_threads >= n || n > AMBIENT_THREADS_LIMIT) return;
        for(size_t i = this->num_threads; i < n; i++){
            pthread_create(&this->pool[i], NULL, &controller::stream, &this->tasks[i]);
        }
        this->num_threads = n;
    }

    inline void controller::allocate_threads(){
        for(size_t i = 1; i < AMBIENT_THREADS_LIMIT; i++){
            this->tasks[i].id = i;
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
                //printf("WARNING: MASTER HAS NULL INSTRUCTIONS! %d\n", (int)this->workload);
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

    inline void controller::alloc_block(revision& r, size_t x, size_t y){
        r.embed(r.spec->alloc(), x, y, r.spec->get_bound());
    }

    inline void controller::calloc_block(revision& r, size_t x, size_t y){
        r.embed(r.spec->calloc(), x, y, r.spec->get_bound());
    }

    // note: ufetch_block is used only by pt_fetch in user-space
    inline revision::entry& controller::ufetch_block(revision& r, size_t x, size_t y){
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

    inline revision::entry& controller::ifetch_block(revision& r, size_t x, size_t y){
        pthread_mutex_lock(&r.mutex);
        if(r.get_generator() == NULL) this->atomic_receive(r, x, y); // check the stacked operations for the block
        pthread_mutex_unlock(&r.mutex);
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
        printf("Warning: using slow block-locking!\n");
        bool acquired = r.block(x,y).trylock();
        return acquired;
    }

    inline void controller::unlock_block(revision& r, size_t x, size_t y){
        r.block(x,y).unlock();
    }

    inline void controller::atomic_receive(revision& r, size_t x, size_t y){
        std::list<cfunctor*>& list = r.block(x,y).get_assignments();
        std::list<cfunctor*>::iterator it = list.begin(); 
        while(it != list.end())
        {
            bool ready = true;
            std::list<revision*>& deps = (*it)->get_dependencies();
            for(std::list<revision*>::iterator d = deps.begin(); d != deps.end(); ++d){
                if(*d == &r) continue;
                pthread_mutex_lock(&(*d)->mutex);
                if((*d)->get_generator() != NULL){
                    (*d)->content.assignments.push_back(*it);
                    ready = false;
                }
                pthread_mutex_unlock(&(*d)->mutex);
                if(!ready) break;
            }
            if(ready) this->execute_mod(*it, dim2(x,y));
            list.erase(it++);
        }
    }

    inline void controller::atomic_complete(cfunctor* op){
        pthread_mutex_lock(&this->mutex);
        this->workload--;
        pthread_mutex_unlock(&this->mutex);

        std::list<revision*>& list = op->get_derivatives();
        for(std::list<revision*>::iterator it = list.begin(); it != list.end(); ++it){
            pthread_mutex_lock(&(*it)->mutex);
            (*it)->reset_generator();
            ambient::controller.atomic_receive(**it, 0, 0);
            pthread_mutex_unlock(&(*it)->mutex);
        }
        delete op;
    }

    inline void controller::mute(){
        this->muted = true;
    }

    inline void controller::unmute(){
        this->muted = false;
    }

    inline void controller::conditional_flush(){
        // you can call flush any time some condition
        // has been satisfied (i.e. memory has ended)
    }

    inline void controller::flush(){
        if(this->stack.empty()) return;
        //printf("PLAYOUT WITH %d\n", (int)this->stack.size());
        while(!this->stack.end_reached())
            this->stack.pick()->logistics();
        this->master_stream(this->tasks);
        this->stack.reset();
    }

/*  DEBUG VERSION:
    inline void controller::flush(){
        static __a_timer ts("ambient: time of small sets");
        static __a_timer tl("ambient: time of large sets");
        __a_timer& time = (this->stack.size() > 100)? ts : tl;
        if(this->stack.empty()) return;

        //while(!this->stack.end_reached()) 
        //    this->stack.pick()->weight();
        //this->stack.sort();

        time.begin();
        while(!this->stack.end_reached())
            this->stack.pick()->logistics();
        this->master_stream(this->tasks);
        this->stack.reset();
        time.end();
    } */
    
    inline packet* package(revision& r, const char* state, int x, int y, int dest){
        void* header = r.block(x,y).get_memory();
        //if(header == NULL) printf("HEADER IS NULL (SWAPPED)\n");
        packet* package = pack(*(packet_t*)r.spec->get_packet_t(), 
                               header, dest, "P2P", r.sid, state, x, y, NULL);
        return package;
    }

    inline void forward_block(packet& cmd){
        /*packet& c = static_cast<packet&>(cmd);
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
        }*/
    }

    inline void accept_block(packet& cmd){
        /*packet& c = static_cast<packet&>(cmd);
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

        ambient::controller.atomic_receive(l, x, y); // calling controller event handlers*/
    }

} } }
