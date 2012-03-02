#ifndef AMBIENT_TASKLIST_H
#define AMBIENT_TASKLIST_H

#include <pthread.h>

#define STACK_CONTENT_RESERVATION 10
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
