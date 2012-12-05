#ifndef AMBIENT_CONTROLLERS_VELVET_CONTROLLER
#define AMBIENT_CONTROLLERS_VELVET_CONTROLLER

#include "ambient/controllers/velvet/cfunctor.h"
#include "ambient/controllers/velvet/context.h"
#include "ambient/controllers/velvet/iteratable.h"
#include "ambient/utils/collector.h"

#include <cilk/cilk.h>

namespace ambient { namespace controllers { namespace velvet {

    using ambient::models::velvet::history;
    using ambient::models::velvet::revision;
    using ambient::models::velvet::memspec;

    class controller : public singleton< controller >
    {
    public:
        controller();
       ~controller();
        void flush();
        void   acquire(channels::mpi::channel* channel);
        void   schedule(cfunctor* op);

        void calloc(revision& r);
        void alloc (revision& r);
        void free  (revision& r);
        revision& ufetch(revision& r);
        void ifetch(revision& r);
        void unlock_revision(revision* arg);
        void unlink_revision(revision* arg);

        template<typename T> void destroy(T* o);

        void atomic_receive(revision& r);
    public:
        std::vector< cfunctor* > chains;
        std::vector< cfunctor* > mirror;
        int arity;
        collector garbage;
    };
    
} } }

namespace ambient {
    extern controllers::velvet::controller& controller;
}

#include "ambient/controllers/velvet/controller.hpp"
#include "ambient/utils/collector.hpp"
#include "ambient/controllers/velvet/context.hpp"
#include "ambient/controllers/velvet/iteratable.hpp"
#include "ambient/controllers/velvet/cfunctor.hpp"
#endif
