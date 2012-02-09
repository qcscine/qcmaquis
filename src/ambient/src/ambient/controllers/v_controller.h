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
        void   acquire(channels::ichannel* channel);
        void   push(models::imodel::modifier* op);

        models::imodel::layout::entry& fetch_block(models::imodel::revision& r, size_t i, size_t j);
        models::imodel::layout::entry& ifetch_block(models::imodel::revision& r, size_t i, size_t j);
        void unlock_revision(models::imodel::revision* arg);
        void unlink_revision(models::imodel::revision* arg);

        void flush();
        void init_threads();
       ~v_controller();

    private:
        models::imodel* model;
        channels::ichannel* channel;
        touchstack< models::imodel::modifier* > stack;
        pthread_t* pool;
        tasklist* tasks;
        dim2 idx;
        workgroup_context context;
    };

} }
#endif
