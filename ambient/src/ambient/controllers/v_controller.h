#ifndef AMBIENT_V_CONTROLLER_H
#define AMBIENT_V_CONTROLLER_H
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
        void   push_mod(mod*);

        models::imodel::layout::entry& fetch_block(models::imodel::revision& r, size_t i, size_t j);
        models::imodel::layout::entry& ifetch_block(models::imodel::revision& r, size_t i, size_t j);
        models::imodel::layout::entry& init_block(models::imodel::revision& r, size_t i, size_t j);
        void unlock_revision(models::imodel::revision* arg);
        void unlink_revision(models::imodel::revision* arg);

        void flush();
        void init_threads();
        void atomic_complete();
        void atomic_receive(models::imodel::revision& r, size_t i, size_t j);
       ~v_controller();

    private:
        pthread_mutex_t mutex;
        models::imodel* model;
        channels::ichannel* channel;
        touchstack< models::imodel::modifier* > stack;
        pthread_t* pool;
        tasklist* tasks;
        size_t workload;
    };

} }
#endif
