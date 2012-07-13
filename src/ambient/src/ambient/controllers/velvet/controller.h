#ifndef AMBIENT_CONTROLLERS_VELVET_CONTROLLER
#define AMBIENT_CONTROLLERS_VELVET_CONTROLLER
#include "ambient/utils/touchstack.h"
#include "ambient/controllers/velvet/cfunctor.h"
#include "ambient/controllers/velvet/context.h"
#include "ambient/controllers/velvet/iteratable.h"
#include "ambient/utils/tasklist.hpp"

namespace ambient { namespace controllers { namespace velvet {

    using ambient::models::velvet::revision;
    using ambient::models::velvet::memspec;

    class controller : public singleton< controller >
    {
    public:
        controller();
        static void* stream(void* list);
        inline void   master_stream(void* list); // specialized version for the main thread
        inline void   acquire(channels::mpi::channel* channel);
        inline void   push(cfunctor* op);
        inline void   execute_mod(cfunctor* op, dim2 pin);

        inline void alloc_block (revision& r, size_t x, size_t y);
        inline void calloc_block(revision& r, size_t x, size_t y);
        inline revision::entry& ufetch_block(revision& r, size_t x, size_t y);
        inline revision::entry& ifetch_block(revision& r, size_t x, size_t y);
        inline void unlock_revision(revision* arg);
        inline void unlink_revision(revision* arg);

        inline void flush();
        inline void conditional_flush();
        inline void allocate_threads();
        inline void set_num_threads(size_t n);
        inline size_t get_num_threads() const;
        inline void atomic_complete(cfunctor* op);
        inline void atomic_receive(revision& r, size_t x, size_t y);
        inline ~controller();
    public:
        bool muted;
    private:
        touchstack< cfunctor* > stack;
        pthread_t pool[AMBIENT_THREADS_LIMIT];
        tasklist_async tasks[AMBIENT_THREADS_LIMIT];
        tasklist_async resolutionq[AMBIENT_THREADS_LIMIT];
        size_t workload;
        size_t num_threads;
        size_t rrn;
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
