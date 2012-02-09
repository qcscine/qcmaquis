#include "ambient/ambient.h"
#include "ambient/controllers/context.h"

extern pthread_key_t pthread_env;

namespace ambient { namespace controllers {

    context::context()
    :grp(NULL)
    { 
    }

    context& context::operator()(const int rank){
        return *this;
    }

    void context::set_group(channels::group* grp){
        this->op->set_group(grp);
        this->grp = grp;
    }

    channels::group* context::get_group(){
        if(this->grp == NULL){
            printf("Attempting to access NULL scope, check if select() was called\n");
        }
        return this->grp;
    }

    int context::get_master_g(){
        return this->grp->get_master_g();
    }

    int context::get_rank(){
        return this->grp->get_rank();
    }

    int context::get_size(){
        return this->grp->get_size();
    }

    dim2 context::get_block_id(){
        void* id = pthread_getspecific(pthread_env);
        assert(id != NULL); 
        return *(dim2*)id;
    }

    void context::set_block_id(dim2 value){
        void* id = pthread_getspecific(pthread_env);
        if(id == NULL){
            id = malloc(sizeof(dim2));
            pthread_setspecific(pthread_env, id);
        }
        *(dim2*)id = value;
    }

    const char* context::get_name(){
        return this->grp->get_name();
    }

    void context::set_op(models::imodel::modifier* op){
        this->op = op;
    }

    models::imodel::modifier* context::get_op(){
        return this->op;
    }

    bool context::involved(){
        if(this->state == MARKUP) return false;
        return this->grp->involved();
    }

    bool context::is_master(){
        return this->grp->is_master();
    }

} }
