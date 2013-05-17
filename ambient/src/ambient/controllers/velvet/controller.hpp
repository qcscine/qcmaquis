#include "ambient/utils/io.hpp"
#include "ambient/utils/timings.hpp"

namespace ambient { namespace controllers { namespace velvet {

    inline controller::~controller(){ }

    inline controller::controller()
    {
        this->stack_m.reserve(AMBIENT_STACK_RESERVE);
        this->stack_s.reserve(AMBIENT_STACK_RESERVE);
        this->chains = &this->stack_m;
        this->mirror = &this->stack_s;
        ambient::channel.init();
        if(ambient::rank()) ambient::rank.mute();
        this->context = new ambient::scope<base>();
    }

    inline bool controller::tunable(){
        return context->tunable();
    }

    template<complexity O>
    inline void controller::schedule(){
        context->toss();
    }

    inline void controller::intend_fetch(history* o){
        revision* r = o->back();
        if(r == NULL) return context->consider_allocation(o->spec.size); // do we have to?
        context->consider_transfer(r->extent, r->state);
    }

    inline void controller::intend_write(history* o){
        context->consider_allocation(o->spec.size);
    }

    inline void controller::set_context(scope* s){
        this->context_c = this->context;
        this->context = s; // no nesting
    }

    inline void controller::pop_context(){
        this->context = this->context_c;
    }

    inline bool controller::remote(){
        return (this->context->state == ambient::remote);
    }

    inline bool controller::local(){
        return (this->context->state == ambient::local);
    }

    inline void controller::flush(){
        typedef typename std::vector<cfunctor*>::const_iterator veci;
        ambient::model.clock++;
        #ifdef AMBIENT_COMPUTATIONAL_DATAFLOW
        printf("ambient::parallel graph dim: %d\n", chains->size());
        for(veci i = chains->begin(); i != chains->end(); ++i)
            printf("op%d[label=\"%s\"]\n", (*i)->id(), (*i)->name());
        #endif
        AMBIENT_SMP_ENABLE
        while(!chains->empty()){
            for(veci i = chains->begin(); i != chains->end(); ++i){
                if((*i)->ready()){
                    cfunctor* task = *i;
                    AMBIENT_THREAD task->invoke();
                    #ifdef AMBIENT_COMPUTATIONAL_DATAFLOW
                    for(int n = 0; n < task->deps.size(); ++n)
                        printf("op%d[label=\"%s\"]\nop%d -> op%d\n", task->deps[n]->id(), task->deps[n]->name(), 
                                                                     task->id(), task->deps[n]->id());
                    #endif
                    mirror->insert(mirror->end(), (*i)->deps.begin(), (*i)->deps.end());
                }else mirror->push_back(*i);
            }
            chains->clear();
            std::swap(chains,mirror);
        }
        AMBIENT_SMP_DISABLE
    }

    inline bool controller::empty(){
        return this->chains->empty();
    }

    inline void controller::clear(){
        this->garbage.clear();
    }

    inline bool controller::queue(cfunctor* f){
        this->chains->push_back(f);
        return true;
    }

    inline void controller::alloc(revision& r){
        r.embed(ambient::pool.malloc(r.extent, r.region));
    }

    inline void controller::calloc(revision& r){
        alloc(r); memset(r.data, 0, r.extent);
    }

    inline void controller::free(revision& r){
        return ambient::pool.free(r.data, r.extent, r.region);
    }

    inline void controller::sync(revision* r){
        if(context->round == 1) return;
        if(ambient::model.common(r)) return;
        if(ambient::model.feeds(r)) ambient::controllers::velvet::set<revision, AMBIENT_NUM_PROCS>::spawn(*r);
        else ambient::controllers::velvet::get<revision, AMBIENT_NUM_PROCS>::spawn(*r);
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
