#include "ambient/utils/io.hpp"
#include "ambient/utils/timings.hpp"

#define CONTROLLER_CHAINS_RESERVE 65536
#define BOUND 0

namespace ambient { namespace controllers { namespace velvet {

    using ambient::channels::mpi::packet_t;
    using ambient::channels::mpi::packet;

    inline controller::~controller(){ }

    inline controller::controller() : depth(0), uniform(false)
    {
        this->stack_m.reserve(CONTROLLER_CHAINS_RESERVE);
        this->stack_s.reserve(CONTROLLER_CHAINS_RESERVE);
        this->chains = &this->stack_m;
        this->mirror = &this->stack_s;
        ambient::channel.init();
        if(ambient::rank()) ambient::rank.mute();
        this->context = new ambient::scope<base>();
    }

    template<complexity O>
    inline void controller::schedule(){
        /*if(uniform) set_state(ambient::common); 
        else if(tuning.decide() == ambient::rank()) set_state(ambient::feed);
        else set_state(ambient::stub);
        tuning.repeat();*/
    }

    /*inline void controller::score_r(history* o){
        revision* r = o->back();
        if(r == NULL) tuning.add_as_new(o->spec.size);
        else if(r->state == ambient::stub) tuning.add_as_remote(r->extent);
        else if(r->state == ambient::feed) tuning.add_as_local(r->extent);
    }

    inline void controller::score_w(history* o){
        tuning.add_as_new(o->spec.size);
    }

    inline void controller::score_rw(history* o){
        score_r(o);
        score_w(o);
    }

    inline void controller::set_state(ambient::rstate s){
        this->state = s;
    }*/

    inline void controller::set_context(scope* s){
        this->depth++;
        if(depth > 1) printf("NESTED SCOPES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n\n");
        this->context_c = this->context;
        this->context = s;
    }

    inline void controller::pop_context(){
        this->depth--;
        this->context = this->context_c;
    }

    inline bool controller::remote(){
        //return (this->state == ambient::stub); //(this->context->state == scope::REMOTE);
        return (this->context->state == ambient::REMOTE);
    }

    inline bool controller::local(){
        //return (this->state == ambient::feed); //(this->context->state == scope::LOCAL);
        return (this->context->state == ambient::LOCAL);
    }

    inline void controller::flush(){
        typedef typename std::vector<cfunctor*>::const_iterator veci;
        ambient::model.clock++;
        //printf("Cumulative imbalance: %.2f\n", ((float)tuning.load[0]/tuning.load[1]));
        //tuning.clear();
        #ifdef AMBIENT_OMP
        #pragma omp parallel
        {
            #pragma omp single
        #endif
            while(!chains->empty()){
                for(veci i = chains->begin(); i != chains->end(); ++i){
                    if((*i)->ready()){
                        cfunctor* task = *i;
                        #ifdef AMBIENT_OMP
                        #pragma omp task untied
                        #endif
                        AMBIENT_THREAD task->invoke();
                        mirror->insert(mirror->end(), (*i)->deps.begin(), (*i)->deps.end());
                    }else mirror->push_back(*i);
                }
                chains->clear();
                std::swap(chains,mirror);
            }
        #ifdef AMBIENT_OMP
        }
        #endif
    }

    inline bool controller::empty(){
        this->chains->empty();
    }

    inline void controller::clear(){
        this->garbage.clear();
    }

    inline bool controller::queue(cfunctor* f){
        this->chains->push_back(f);
        return true;
    }

    inline void controller::alloc(revision& r){
        r.embed(ambient::pool.malloc(r.extent + BOUND, r.region), BOUND);
    }

    inline void controller::calloc(revision& r){
        alloc(r); memset(r.data, 0, r.extent);
    }

    inline void controller::free(revision& r){
        return ambient::pool.free(r.header, r.extent + BOUND, r.region);
    }

    inline void controller::sync(revision* r){
        if(context->round == 1) return; // for single process version
        if(ambient::model.common(r)) return;
        if(ambient::model.feeds(r)) ambient::controllers::velvet::set<revision>::spawn(*r) >> (1-ambient::rank()); // for 2 processes
        else ambient::controllers::velvet::get<revision>::spawn(*r);
    }

    inline void controller::lsync(revision* r){
        if(ambient::model.common(r)) return;
        if(!ambient::model.feeds(r)) ambient::controllers::velvet::get<revision>::spawn(*r);
    }

    inline void controller::rsync(revision* r){
        if(ambient::model.common(r)) return;
        if(ambient::model.feeds(r)) ambient::controllers::velvet::set<revision>::spawn(*r) >> context->sector;
    }

    inline void controller::lsync(transformable* v){
        if(context->round == 1) return;
        ambient::controllers::velvet::set<transformable, AMBIENT_NUM_PROCS>::spawn(*v);
    }

    inline void controller::rsync(transformable* v){
        ambient::controllers::velvet::get<transformable, AMBIENT_NUM_PROCS>::spawn(*v);
    }

    template<typename T> void controller::destroy(T* o){
        this->garbage.push_back(o);
    }

    inline void controller::persist(history* o){
        o->back()->region = PERSIST_REGION;
    }

} } }
