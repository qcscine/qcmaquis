#include "ambient/utils/tasklist.h"

namespace ambient{

    tasklist::tasklist()
    : seed(NULL), tail(NULL), active(true)
    {
        pthread_mutex_init(&this->mutex, NULL);
    }
    
    tasklist::~tasklist(){
        pthread_mutex_destroy(&this->mutex);
    }

    void tasklist::add_task(void* content){ // delegate thread
        pthread_mutex_lock(&this->mutex);

        if(this->tail == NULL){
            this->tail = new task(content);
            this->seed = this->tail;
        }else{
            this->tail->set_next(new task(content));
            this->tail = this->tail->next();
        }

        pthread_mutex_unlock(&this->mutex);
    }

    void* tasklist::get_task(){ // worker thread
        if(this->seed == NULL) return NULL;
        else if(this->seed == this->tail){
            pthread_mutex_lock(&this->mutex);
            void* data = this->seed->content;
            delete this->seed;
            this->seed = NULL;
            this->tail = NULL;
            pthread_mutex_unlock(&this->mutex);
            return data;
        }
        task* seed = this->seed;
        void* data = seed->content;
        this->seed = seed->next();
        delete seed;
        return data;
    }

    tasklist::task::task(void* data)
    : content(data), pair(NULL)
    {
    }

    void tasklist::task::set_next(task* t){
        this->pair = t;
    }

    tasklist::task* tasklist::task::next(){
        return this->pair;
    }


} // namespace ambient
