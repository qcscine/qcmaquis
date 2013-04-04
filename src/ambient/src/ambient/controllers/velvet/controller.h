#ifndef AMBIENT_CONTROLLERS_VELVET_CONTROLLER
#define AMBIENT_CONTROLLERS_VELVET_CONTROLLER

#include "ambient/controllers/velvet/cfunctor.h"
#include "ambient/controllers/velvet/iteratable.h"
#include "ambient/utils/collector.h"

namespace ambient { namespace controllers { namespace velvet {

    using ambient::models::velvet::revision;

    class controller : public singleton< controller >
    {
    public:
        class scope {
        public:
            int sector;
            int round;
            int gauge;
            ambient::locality state;
            virtual bool tunable() = 0;
            virtual void consider_transfer(size_t size, ambient::locality l){}
            virtual void consider_allocation(size_t size){}
            virtual void toss(){}
        };

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

        bool tunable();
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
