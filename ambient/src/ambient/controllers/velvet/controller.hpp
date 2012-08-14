#include "ambient/utils/io.hpp"
#include "ambient/utils/timings.hpp"

namespace ambient { namespace controllers { namespace velvet {

    using ambient::channels::mpi::packet_t;
    using ambient::channels::mpi::packet;

    inline controller::~controller(){
        delete[] this->stacks;
    }

    inline controller::controller()
    : muted(false)
    {
        this->acquire(&ambient::channel);
        this->stacks = new touchstack< cfunctor* >[__cilkrts_get_nworkers()];
        this->garbage = new touchstack< history* >[__cilkrts_get_nworkers()];
        this->mgarbage = new touchstack< void* >[__cilkrts_get_nworkers()];
    }
        
    inline void controller::destroy(history* o){
        this->garbage[__cilkrts_get_worker_number()].push_back(o);
    }
        
    inline void controller::deallocate(void* o){
        this->mgarbage[__cilkrts_get_worker_number()].push_back(o);
    }

    inline void controller::push(cfunctor* op){
        this->stacks[__cilkrts_get_worker_number()].push_back(op);
    }

    inline void controller::execute(chain* op){
        op->execute();
        delete op;
    }

    inline void controller::flush(){
        for(int i = 0; i < __cilkrts_get_nworkers(); i++){
            while(!this->stacks[i].end_reached())
                this->stacks[i].pick()->logistics();
            this->stacks[i].reset();
        }
        if(this->chains.size() == 1){
            this->execute((chain*)this->chains.pick());
            this->chains.reset();
        }else{
            touchstack< chain* >* chains = &this->chains;
            touchstack< chain* >* mirror = &this->mirror;

            while(!chains->empty()){
                while(!chains->end_reached()){
                    chain* op = chains->pick();
                    if(op->ready()) cilk_spawn this->execute(op);
                    else mirror->push_back(op);
                }
                chains->reset();
                std::swap(chains,mirror);
            }
            cilk_sync;
        }
        ambient::bulk_pool.refresh();
        for(int i = 0; i < __cilkrts_get_nworkers(); i++){
            garbage[i].clear();
            mgarbage[i].purge();
        }
        // implicit cilk_sync //
    }

    inline void controller::execute_mod(cfunctor* op){
        if(!this->chains.empty() && op->match(this->chains.back())) 
            this->chains.back()->push_back(op);
        else
            this->chains.push_back(new chain(op));
    }

    inline void controller::atomic_receive(revision& r){
        /*std::list<cfunctor*>& list = r.content.assignments;
        std::list<cfunctor*>::iterator it = list.begin(); 
        while(it != list.end()){
            this->execute_mod(*it);
            list.erase(it++);
        }*/ // should be rewritten for MPI
    }

    inline void controller::acquire(channels::mpi::channel* channel){
        channel->init();
    }

    inline void controller::alloc(revision& r){
        r.embed(r.spec->alloc(), r.spec->get_bound());
    }

    inline void controller::calloc(revision& r){
        r.embed(r.spec->calloc(), r.spec->get_bound());
    }

    inline revision& controller::ufetch(revision& r){
        return r;
    }

    inline void controller::ifetch(revision& r){
        //this->atomic_receive(r);
    }

    inline void controller::conditional_flush(){
        // you can call flush any time some condition
        // has been satisfied (i.e. memory has ended)
    }

    inline void forward(packet& cmd){ }

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
    
    /* MPI required parts:
    inline packet* package(revision& r, const char* state, int x, int y, int dest){
        void* header = r.block(x,y).get_memory();
        //if(header == NULL) printf("HEADER IS NULL (SWAPPED)\n");
        packet* package = pack(*(packet_t*)r.spec->get_packet_t(), 
                               header, dest, "P2P", r.sid, state, x, y, NULL);
        return package;
    }

    inline void forward(packet& cmd){
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
            forward(cmd);             // and forwarding
        }else{
            l.get(x,y).get_path().push_back(c.get<int>(A_LAYOUT_P_OWNER_FIELD));
        }
    }

    inline void accept(packet& cmd){
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
    } */

    /* Full fledged heavy version:
    inline revision::entry& controller::ifetch(revision& r){
        assert(r.get_placement() != NULL);
        if(r.block(x,y).valid())
            this->atomic_receive(r.get_layout(), x, y);
        else if(r.get_placement()->is_master()){
            if(r.get_layout().marked(x, y)){
                this->alloc_block(r.get_layout(), x, y); // ! making all allocations inside computational kernels
                this->atomic_receive(r.get_layout(), x, y);
            }else{
                printf("UNEXPECTED ERROR IN IFETCH BLOCK -- asking for %lu %lu!\n", x, y);
                if(r.get_generator()->get_group() != NULL) // we already know the generation place
                    ambient::channel.ifetch(r.get_generator()->get_group(), r.get_layout().id(), x, y);
             // else
             //     leaving on the side of the generator to deliver (we don't really know who generates but the destination is already known)
            }
        }else if(!r.block(x,y).valid() && !r.block(x,y).requested()){
            ambient::channel.ifetch(r.get_placement(), r.get_layout().id(), x, y);
        }
        return r.block(x,y);
    } 

    // Full fledged heavy version: (note that ufetch_block is used only by pt_fetch in user-space)
    inline revision::entry& controller::ufetch_block(revision& r, size_t x, size_t y){
        if(r.block(x,y).valid()){
            return r.block(x,y);
        }else if(r.get_placement() == NULL || r.get_placement()->is_master()){
            assert(r.get_generator()->get_group() != NULL);
            if(r.get_generator()->get_group()->is_master()){
              return this->alloc_block(r.get_layout(), x, y);
            }
        }else if(!r.block(x,y).requested()){
            ambient::channel.ifetch(r.get_placement(), *r.get_layout().id(), x, y);
        } // blocks should be already requested
        return r.block(x,y);
    } */

} } }
