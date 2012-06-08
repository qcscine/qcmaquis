#ifndef AMBIENT_UTILS_TASKLIST
#define AMBIENT_UTILS_TASKLIST
#define STACK_CONTENT_RESERVATION 10

namespace ambient{

    class tasklist {
    public:
        class task {
        public:
            inline task(void* data):content(data),pair(NULL){}
            inline void set_next(task* t){ pair = t; }
            inline task* next(){ return pair; }
            void* content;
            task* pair;
        };

        inline tasklist()
        : seed(NULL), tail(NULL), active(true), idle(false)
        {
            pthread_mutex_init(&mutex, NULL);
        }
        inline ~tasklist(){
            pthread_mutex_destroy(&mutex);
        }
        inline void add_task(void* content){ // delegate thread
            pthread_mutex_lock(&mutex);
        
            if(this->seed == NULL){
                this->seed = new task(content);
                this->tail = this->seed;
            }else{
                this->tail->set_next(new task(content));
                this->tail = this->tail->next();
            }
        
            pthread_mutex_unlock(&mutex);
        }
        inline void* get_task(){ // worker thread
            void* data = NULL;
            pthread_mutex_lock(&mutex);
            if(this->seed != NULL){
                data = this->seed->content;
                if(this->seed == this->tail){
                    delete this->seed;
                    this->seed = NULL;
                    this->tail = NULL;
                }else{
                    task* seed = this->seed;
                    this->seed = this->seed->next();
                    delete seed;
                }
            }
            pthread_mutex_unlock(&mutex);
            return data;
        }

        task* seed;
        task* tail;
        pthread_mutex_t mutex;
        bool active;
        bool idle;
        size_t id;
    };

}

#endif
