#include "ambient/ambient.h"
#include "ambient/controllers/context.h"

extern pthread_key_t pthread_env;
extern pthread_key_t pthread_tid;

namespace ambient { namespace controllers {

    context::context()
    :grp(NULL), state(MARKUP)
    { 
    }

    void context::set_group(channels::group* grp){
        this->op->set_group(grp);
        this->grp = grp;
    }

    void context::set_tid(size_t value){
        void* tid = pthread_getspecific(pthread_tid);
        if(tid == NULL){
            tid = malloc(sizeof(size_t));
            pthread_setspecific(pthread_tid, tid);
        }
        *(size_t*)tid = value;
    }

    void context::set_block_id(dim2 value){
        void* id = pthread_getspecific(pthread_env);
        if(id == NULL){
            id = malloc(sizeof(dim2));
            pthread_setspecific(pthread_env, id);
        }
        *(dim2*)id = value;
    }

    size_t context::get_revision_base(const models::imodel::object* o){
        if(this->state == MARKUP) return o->get_revision_base();
        return o->get_thread_revision_base();
    }

    void context::set_revision_base(models::imodel::object* o, size_t base){
        if(this->state == MARKUP) o->set_revision_base(base);
        else o->set_thread_revision_base(base);
    }

} }
