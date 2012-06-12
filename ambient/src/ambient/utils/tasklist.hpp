#ifndef AMBIENT_UTILS_TASKLIST
#define AMBIENT_UTILS_TASKLIST
#define STACK_CONTENT_RESERVATION 10

namespace ambient{

    using controllers::velvet::cfunctor;

    class tasklist {
    public:
        class task {
        public:
            inline task(cfunctor* f, dim2 pin):f(f),pin(pin),pair(NULL){}
            inline task* next(task* t){ pair = t; return t; }
            inline task* next(){ return pair; }
            task* pair;
            cfunctor* f;
            dim2 pin;
        };

        inline tasklist()
        : seed(NULL), tail(NULL), active(true)
        {
            pthread_mutex_init(&mutex, NULL);
        }
        inline ~tasklist(){
            pthread_mutex_destroy(&mutex);
        }
        inline void add_task(task* t){ // delegate thread
            pthread_mutex_lock(&mutex);
            if(this->seed == NULL) this->tail = this->seed = t;
            else this->tail = this->tail->next(t);
            pthread_mutex_unlock(&mutex);
        }
        inline task* get_task(){ // worker thread
            pthread_mutex_lock(&mutex);
            task* data = this->seed;
            if(this->seed != NULL){
                if(this->seed == this->tail) this->tail = this->seed = NULL;
                else this->seed = this->seed->next();
            }
            pthread_mutex_unlock(&mutex);
            return data;
        }

        task* seed;
        task* tail;
        pthread_mutex_t mutex;
        bool active;
        size_t id;
    };

}

#endif
