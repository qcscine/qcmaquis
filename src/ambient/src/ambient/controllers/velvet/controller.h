#ifndef AMBIENT_CONTROLLERS_VELVET_CONTROLLER
#define AMBIENT_CONTROLLERS_VELVET_CONTROLLER

#include "ambient/utils/touchstack.h"
#include "ambient/controllers/velvet/cfunctor.h"
#include "ambient/controllers/velvet/chain.h"
#include "ambient/controllers/velvet/context.h"
#include "ambient/controllers/velvet/iteratable.h"
#include "ambient/utils/tasklist.hpp"

#include <cilk/cilk.h>

namespace ambient { namespace controllers { namespace velvet {

    using ambient::models::velvet::history;
    using ambient::models::velvet::revision;
    using ambient::models::velvet::memspec;

    class controller : public singleton< controller >
    {
    public:
        controller();
        inline void   acquire(channels::mpi::channel* channel);
        inline void   push(cfunctor* op);
        inline void   execute_mod(cfunctor* op);

        inline void alloc (revision& r);
        inline void calloc(revision& r);
        inline revision& ufetch(revision& r);
        inline void ifetch(revision& r);
        inline void unlock_revision(revision* arg);
        inline void unlink_revision(revision* arg);
        inline void destroy(history* o);
        inline void deallocate(void* o);

        inline void flush();
        inline void conditional_flush();
        inline void atomic_receive(revision& r);
        inline ~controller();
    public:
        bool muted;
    private:
        touchstack< history* >* garbage;
        touchstack< void* >* mgarbage;
        touchstack< cfunctor* >* stacks;
        inline void execute(chain* op);
        touchstack< chain* > chains;
        touchstack< chain* > mirror;
    };
    
} } }

namespace ambient {
    extern controllers::velvet::controller& controller;
}

#include "ambient/controllers/velvet/controller.hpp"
#include "ambient/controllers/velvet/context.hpp"
#include "ambient/controllers/velvet/iteratable.hpp"
#include "ambient/controllers/velvet/chain.hpp"
#include "ambient/controllers/velvet/cfunctor.hpp"
#endif
