#include "ambient/ambient.h"
#include "ambient/model.h"

namespace ambient{
// model //
    data_model& data_model::instance(){
        static data_model* singleton = NULL;
        if(!singleton) singleton = new data_model();
        return *singleton;
    }

    data_model::data_model(){

    }

    p_object* data_model::get_object(unsigned int* hash, unsigned int hash_len, unsigned int id) const {
        return this->map.find(hash, hash_len, id)->object;
    }

    void data_model::add_object(unsigned int* hash, unsigned int hash_len, p_object* object){
        object->layout = new core::layout_table(object);
        object->group_id = hash;
        object->id = this->map.insert(hash, hash_len, object->layout);
    }

    void data_model::update_object(groups::group* scope, p_object* object){
        if(object->state == SERIAL) return;
        if(!object->id) this->add_object(scope->id.first, scope->id.second, object);
        if(object->state == ABSTRACT)       object->state = COMPOSING;
        else if(object->state == COMPOSING) object->state = GENERIC;
        if(!object->consted){
            object->timestamp++;
            object->set_scope(scope);
        }
    }
// model //

// controller //
    d_controller& d_controller::instance(){
        static d_controller* singleton = NULL;
        if(!singleton) singleton = new d_controller();
        return *singleton;
    }

    d_controller::d_controller(){

    }

    void d_controller::acquire(i_channel* channel){
        this->channel = channel;
    }
// controller //

// channel // 
    mpi_channel& mpi_channel::instance(){
        static mpi_channel* singleton = NULL;
        if(!singleton) singleton = new mpi_channel();
        return *singleton;
    }

    mpi_channel::mpi_channel(){
        this->controller = &ambient::controller;
        this->controller->acquire(this);
        
    }
// channel // 
}
