#ifndef AMBIENT_UTILS_TASKLIST
#define AMBIENT_UTILS_TASKLIST
#define STACK_CONTENT_RESERVATION 10
#include <pthread.h>

namespace ambient{

    class tasklist {
    public:
        class task {
        public:
            task(void* data);
            void set_next(task* t);
            task* next();
            void* content;
            task* pair;
        };

        tasklist();
       ~tasklist();
        void add_task(void* content);
        void* get_task();

        task* seed;
        task* tail;
        pthread_mutex_t mutex;
        bool active;
        bool idle;
        size_t id;
    };

}

#endif
