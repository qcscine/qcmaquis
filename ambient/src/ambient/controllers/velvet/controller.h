#ifndef AMBIENT_CONTROLLERS_VELVET_CONTROLLER
#define AMBIENT_CONTROLLERS_VELVET_CONTROLLER

#include "ambient/controllers/velvet/cfunctor.h"
#include "ambient/controllers/velvet/iteratable.h"
#include "ambient/utils/collector.h"
#include <stack>

namespace ambient { namespace controllers { namespace velvet {

    using ambient::models::velvet::revision;

    class controller : public singleton< controller >
    {
    public:
        class scope {
        public:
            ambient::locality state;
            int sector;
            int round;
            int gauge;
        };

        /*class info {
        public:
            info(){
                footprint = remote = local =
                pin = load[0] = load[1] = rank = 0;
            }
            void repeat(){
                pin = remote = local =
                footprint = 0;
            }
            void clear(){
                load[0] = load[1] = 0;
            }
            int decide(){
                if(power != complexity::N3 && local != remote){
                    if(local > remote){ load[ambient::rank()] += footprint; return ambient::rank(); }
                    else{ load[1-ambient::rank()] += footprint; return (1-ambient::rank()); }
                }

                if(load[0] / std::max(1,load[1]) > 2) rank = 1;
                if(load[1] / std::max(1,load[0]) > 2) rank = 0;
                load[rank] += footprint; //std::pow(pin, power);

                return rank;
            }
            void add_as_new(size_t size){
                footprint += size;
                pin = std::max(pin, size);
            }
            void add_as_local(size_t size){
                pin = std::max(pin, size);
                local += size;
            }
            void add_as_remote(size_t size){
                pin = std::max(pin, size);
                remote += size;
            }
            size_t footprint;
            size_t remote;
            size_t local;
            size_t pin;
            size_t power;
            int load[2];
            int rank;
        } tuning;*/

        controller();
       ~controller();
        bool empty();
        void flush();
        void clear();
        bool queue (cfunctor* f);
        void calloc(revision& r);
        void alloc (revision& r);
        void free  (revision& r);
        void sync  (revision* r);
        void lsync (revision* r);
        void rsync (revision* r);
        void lsync (transformable* v);
        void rsync (transformable* v);
        template<typename T> void destroy(T* o);
        void persist(history* o);

        template<complexity O> void schedule();
        void intend_fetch(history* o);
        void intend_write(history* o);

        void set_context(scope* s);
        void pop_context();
        bool remote();
        bool local();

        scope* context_c;
        scope* context;
    private:
        std::vector< cfunctor* > stack_m;
        std::vector< cfunctor* > stack_s;
        std::vector< cfunctor* >* chains;
        std::vector< cfunctor* >* mirror;
        ambient::collector garbage;
    };
    
} } }

namespace ambient {
    extern controllers::velvet::controller& controller;
}

#include "ambient/controllers/velvet/scope.hpp"
#include "ambient/controllers/velvet/controller.hpp"
#include "ambient/utils/collector.hpp"
#include "ambient/controllers/velvet/iteratable.hpp"
#include "ambient/controllers/velvet/cfunctor.hpp"
#endif
