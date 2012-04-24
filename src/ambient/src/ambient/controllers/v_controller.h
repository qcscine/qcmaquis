#ifndef AMBIENT_CONTROLLERS_V_CONTROLLER
#define AMBIENT_CONTROLLERS_V_CONTROLLER
#include "ambient/controllers/icontroller.h"
#include "ambient/utils/singleton.hpp"
#include "ambient/utils/touchstack.h"
#include "ambient/utils/tasklist.h"

namespace ambient { namespace controllers {

    class v_controller : public icontroller, public singleton< v_controller >
    {
    public: 
        class mod {
            public:
            mod(models::imodel::modifier* m, dim2 pin);
            models::imodel::modifier* m;
            dim2 pin;
        };
        v_controller();
        static void* stream(void* list);
        void   master_stream(void* list); // specialized version for the main thread
        void   acquire(channels::ichannel* channel);
        void   push(models::imodel::modifier* op);
        void   execute_mod(models::imodel::modifier* op, dim2 pin);
        void   execute_free_mod(models::imodel::modifier* op);

        models::imodel::layout::entry* alloc_block(models::imodel::revision& r);
        models::imodel::layout::entry& alloc_block(models::imodel::revision& r, size_t i, size_t j);
        models::imodel::layout::entry& ufetch_block(models::imodel::revision& r, size_t i, size_t j);
        models::imodel::layout::entry& ifetch_block(models::imodel::revision& r, size_t i, size_t j);
        models::imodel::layout::entry& init_block(models::imodel::revision& r, size_t i, size_t j);
        void unlock_revision(models::imodel::revision* arg);
        void unlink_revision(models::imodel::revision* arg);

        void flush();
        void allocate_threads();
        void set_num_threads(size_t n);
        size_t get_num_threads() const;
        void atomic_complete();
        void atomic_receive(models::imodel::revision& r, size_t i, size_t j);

        pthread_mutex_t* get_pool_control_mutex();
       ~v_controller();

    private:
        pthread_mutex_t pool_control_mutex;
        pthread_mutex_t mutex;
        models::imodel* model;
        channels::ichannel* channel;
        touchstack< models::imodel::modifier* > stack;
        pthread_mutex_t* mpool;
        pthread_t* pool;
        tasklist* tasks;
        size_t workload;
        size_t num_threads;
        size_t rrn;
    };

} }

#endif
