#ifndef AMBIENT_CONTROLLERS_VELVET_CONTROLLER
#define AMBIENT_CONTROLLERS_VELVET_CONTROLLER

#include "ambient/controllers/velvet/cfunctor.h"
#include "ambient/controllers/velvet/context.h"
#include "ambient/controllers/velvet/iteratable.h"
#include "ambient/utils/collector.hpp"

#include <cilk/cilk.h>

namespace ambient { namespace controllers { namespace velvet {

    using ambient::models::velvet::history;
    using ambient::models::velvet::revision;
    using ambient::models::velvet::memspec;

    class controller : public singleton< controller >
    {
    public:
        controller();
        void   acquire(channels::mpi::channel* channel);
        void   schedule(cfunctor* op);

        void alloc (revision& r);
        void calloc(revision& r);
        revision& ufetch(revision& r);
        void ifetch(revision& r);
        void unlock_revision(revision* arg);
        void unlink_revision(revision* arg);

        template<typename T> void destroy(T* o);

        void flush();
        void atomic_receive(revision& r);
        ~controller();
    public:
        bool muted;
        collector garbage;
    private:
        std::vector< cfunctor* > chains;
        std::vector< cfunctor* > mirror;
        int arity;
    };
    
} } }

namespace ambient {
    extern controllers::velvet::controller& controller;
}

#include "ambient/controllers/velvet/controller.hpp"
#include "ambient/controllers/velvet/context.hpp"
#include "ambient/controllers/velvet/iteratable.hpp"
#include "ambient/controllers/velvet/cfunctor.hpp"
#endif
