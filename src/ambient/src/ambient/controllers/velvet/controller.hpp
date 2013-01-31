#include "ambient/utils/io.hpp"
#include "ambient/utils/timings.hpp"

#define CONTROLLER_CHAINS_RESERVE 65536
#define BOUND 0

namespace ambient { namespace controllers { namespace velvet {

    using ambient::channels::mpi::packet_t;
    using ambient::channels::mpi::packet;

    inline controller::~controller(){ }

    inline controller::controller(){
        this->stack_m.reserve(CONTROLLER_CHAINS_RESERVE);
        this->stack_s.reserve(CONTROLLER_CHAINS_RESERVE);
        this->chains = &this->stack_m;
        this->mirror = &this->stack_s;
        ambient::channel.init();
        if(ambient::rank() != 0) 
            ambient::rank.verbose = false;
    }

    inline void controller::flush(){
        typedef typename std::vector<cfunctor*>::const_iterator veci;
        ambient::model.clock++;

        while(!chains->empty()){
            for(veci i = chains->begin(); i != chains->end(); ++i){
                if((*i)->ready()){
                    cilk_spawn (*i)->invoke();
                    mirror->insert(mirror->end(), (*i)->deps.begin(), (*i)->deps.end());
                }else mirror->push_back(*i);
            }
            chains->clear();
            std::swap(chains,mirror);
        }
    }

    inline void controller::clear(){
        this->garbage.clear();
    }

    inline void controller::submit(cfunctor* f){
        this->chains->push_back(f);
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

    inline void controller::sync(revision& r, size_t target){
        /*if(ambient::rank() == target && !r.remote()) return; // no-repeat
        if(ambient::rank() != target && !r.valid())  return;
        
        if(r.remote()){ alloc(r); r.state = revision::WAIT; }
        channel.replicate(channels::mpi::vbp::instance(r.header, r.extent, target));*/
    }

    inline void controller::sync(revision& r){
        /*if(!r.origin()){ alloc(r); r.state = revision::WAIT; }
        channel.broadcast( channels::mpi::vbp::instance(r.header, r.extent, ambient::rank()), r.origin() );*/
    }

    template<typename T> void controller::destroy(T* o){
        this->garbage.push_back(o);
    }

} } }
