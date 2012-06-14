#ifndef AMBIENT_UTILS_TASKLIST
#define AMBIENT_UTILS_TASKLIST

namespace ambient{

    using controllers::velvet::cfunctor;

    class tasklist {
    public:
        class task {
        public:
            inline task(cfunctor* f, dim2 pin):f(f),pin(pin),next(NULL){}
            task* next;
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
#ifdef AMBIENT_THREADS
            pthread_mutex_lock(&mutex);
#endif
            if(this->seed == NULL) this->seed = t;
            else this->tail->next = t;
            this->tail = t;
#ifdef AMBIENT_THREADS
            pthread_mutex_unlock(&mutex);
#endif
        }
        inline task* get_task(){ // worker thread
#ifdef AMBIENT_THREADS
            pthread_mutex_lock(&mutex);
#endif
            task* data = this->seed;
            if(this->seed != NULL){
                if(this->seed == this->tail) this->tail = this->seed = NULL;
                else this->seed = this->seed->next;
            }
#ifdef AMBIENT_THREADS
            pthread_mutex_unlock(&mutex);
#endif
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
