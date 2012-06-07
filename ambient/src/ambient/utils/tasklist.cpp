#include "ambient/utils/tasklist.h"

namespace ambient{

    tasklist::tasklist()
    : seed(NULL), tail(NULL), active(true), idle(false)
    {
        pthread_mutex_init(&this->mutex, NULL);
    }
    
    tasklist::~tasklist(){
        pthread_mutex_destroy(&this->mutex);
    }

    void tasklist::add_task(void* content){ // delegate thread
        pthread_mutex_lock(&this->mutex);

        if(this->seed == NULL){
            this->seed = new task(content);
            this->tail = this->seed;
        }else{
            this->tail->set_next(new task(content));
            this->tail = this->tail->next();
        }

        pthread_mutex_unlock(&this->mutex);
    }

    void* tasklist::get_task(){ // worker thread
        void* data = NULL;
        pthread_mutex_lock(&this->mutex);
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
        pthread_mutex_unlock(&this->mutex);
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

}
