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
        *l >> this->mem_dim, this->work_dim, this->item_dim;
        obj->add_revision(l);
    }

    void v_model::update_revision(imodel::revision* r, channels::group* placement){

    }

    v_model::revision* v_model::get_revision(size_t* hash, size_t hash_len, size_t id) const {
        return ((layout*)this->map.find(hash, hash_len, id))->revision;
    }

    v_model& v_model::operator>>(dim2 mem_dim){
        this->mem_dim  = mem_dim;
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

    // {{{ free functions for mangling the data //

    void* solidify(imodel::revision& r){
        imodel::layout& l = r.get_layout();
        size_t iterator = 0;
        char* memory = NULL;
         
        for(size_t j=0; j < l.get_mem_grid_dim().x; j++)
            for(size_t jj=0; jj < l.get_mem_dim().x; jj++)
                for(size_t i=0; i < l.get_mem_grid_dim().y; i++){
                    if(r(i,j).valid()){
                        memory = (char*)realloc(memory, (iterator+1)*l.get_mem_lda());
                        memcpy(memory+iterator*l.get_mem_lda(),                   // copy to
                               &((char*)r(i,j).get_memory())[jj*l.get_mem_lda()], // copy from
                               l.get_mem_lda());                                  // of size
                        iterator++;
                    }
                }
        return memory;
    }

    void disperse(void* data, imodel::revision& r){
        imodel::layout& l = r.get_layout();
        size_t iterator = 0;
        char* memory = (char*)data;

        for(size_t j=0; j < l.get_mem_grid_dim().x; j++)
            for(size_t jj=0; jj < l.get_mem_dim().x; jj++)
                for(size_t i=0; i < l.get_mem_grid_dim().y; i++){
                    if(r(i,j).valid()){
                        memcpy(&((char*)r(i,j).get_memory())[jj*l.get_mem_lda()], // copy from
                               memory+iterator*l.get_mem_lda(),                   // copy to
                               l.get_mem_lda());                                  // of size
                        iterator++;
                    }
                }
        free(data);
    }

    void reduce(v_model::modifier* r){

    }

    void init(v_model::modifier* i){

    }

    // }}}

} }
