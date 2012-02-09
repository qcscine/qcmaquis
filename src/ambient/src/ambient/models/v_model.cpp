#include "ambient/ambient.h"
#include "ambient/models/v_model.h"
#include "ambient/utils/hashmap.h"

namespace ambient { namespace models {

    v_model::v_model()
    : item_dim(dim2(16,16)) // to revert to 128,128
    {
    }

    v_model::~v_model(){
    }

    void v_model::add_revision(imodel::object* obj){
        v_model::layout* l = new v_model::layout(obj->get_dim(), obj->get_t_size());
        *const_cast<size_t*>(&l->sid) = this->map.insert(channel.id().first, channel.id().second, l);
        *const_cast<const size_t **>(&l->gid) = channel.id().first;
        obj->add_revision(l);
    }

    void v_model::update_revision(channels::group* placement, imodel::revision* r){
        //if(r->state == SERIAL) return;
        //if(!r->id) this->add_revision(scope->id.first, scope->id.second, r);
        //if(r->state == ABSTRACT)       r->state = COMPOSING;
        //else if(r->state == COMPOSING) r->state = GENERIC;
        //if(!r->consted){
        //    r->timestamp++;
        //    r->set_scope(scope);
        //}
    }

    v_model::revision* v_model::get_revision(size_t* hash, size_t hash_len, size_t id) const {
        return ((layout*)this->map.find(hash, hash_len, id))->revision;
    }

    v_model& v_model::operator>>(dim2 mem_dim){
        this->mem_dim  = mem_dim;
        //this->default_data_packet_t = new block_packet_t(this->mem_dim*this->item_dim); // to redo in future?
        //this->default_data_packet_t->commit();
        //if(!channel->subscribed(*this->default_data_packet_t)){
        //    channel->subscribe(*this->default_data_packet_t);
        //    channel->add_handler(*this->default_data_packet_t, new core::operation(accept_block, 
        //        channel->get_pipe(*this->default_data_packet_t, channels::ichannel::IN)) );
        //}
        this->work_dim = NULL;
        this->item_dim = NULL;
        return *this;
    }

    v_model& v_model::operator,(dim2 dim){
        if(this->work_dim == NULL)
            this->work_dim = dim;
        else if(this->item_dim == NULL)
            this->item_dim = dim;
        return *this;
    }

    dim2 v_model::get_mem_dim() const {
        return this->mem_dim;  
    }

    dim2 v_model::get_work_dim() const {
        return this->work_dim; 
    }

    dim2 v_model::get_item_dim() const {
        return this->item_dim; 
    }

    // {{{ free functions for mangling the data //

    void* solidify(imodel::revision& instance){ // to correct
        void* data;
        int jumper = 0;
// let's find the solid_lda
        /*
        std::sort(entries.begin(), entries.end(), matrix_order_predicate);
        instance->solid_lda = 0;
        for(int k=0; k < entries.size(); k++) if(entries[k].j == entries[0].j) instance->solid_lda++;
        data = malloc(__a_ceil(std::max(instance->layout->segment_count, instance->layout->request_count)/instance->solid_lda) *
                      instance->get_mem_t_dim().x*instance->get_mem_t_dim().y*instance->t_size*instance->solid_lda);
        assert(data != NULL);
        char* memory = (char*)data;
// let's copy from separate memory blocks to the general one
        for(int k=0; k < entries.size(); k++){
            for(int j=0; j < instance->get_mem_t_dim().x; j++){
                memcpy(memory+j*instance->solid_lda*instance->get_mem_t_dim().y*instance->t_size, 
                       &((char*)instance->block(entries[k].i, entries[k].j)->data)[j*instance->get_block_lda()], 
                       instance->get_mem_t_dim().y*instance->t_size);
            }
            memory += instance->get_mem_t_dim().y*instance->t_size;
            if(++jumper == instance->solid_lda){ jumper = 0; memory += instance->get_mem_t_dim().y*instance->solid_lda*instance->t_size*(instance->get_mem_t_dim().x-1); }
        }*/
        return data;
    }

    void disperse(void* data, imodel::revision& instance){  // to correct
        int jumper = 0;
        /*
        std::sort(entries.begin(), entries.end(), matrix_order_predicate);
        char* memory = (char*)data;
        for(int k=0; k < entries.size(); k++){
            for(int j=0; j < instance->get_mem_t_dim().x; j++){
                 memcpy(&((char*)instance->block(entries[k].i, entries[k].j)->data)[j*instance->get_block_lda()], 
                        memory+j*instance->solid_lda*instance->get_mem_t_dim().y*instance->t_size, 
                        instance->get_mem_t_dim().y*instance->t_size);
            }
            memory += instance->get_mem_t_dim().y*instance->t_size;
            if(++jumper == instance->solid_lda){ jumper = 0; memory += instance->get_mem_t_dim().y*instance->solid_lda*instance->t_size*(instance->get_mem_t_dim().x-1); }
        }
        */
        free(data);
    }

    void reduce(v_model::modifier* r){

    }

    void init(v_model::modifier* i){

    }

    // }}}

} }
