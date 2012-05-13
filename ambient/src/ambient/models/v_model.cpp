#include "ambient/ambient.h"
#include "ambient/models/v_model.h"
#include "ambient/utils/hashmap.h"

namespace ambient { namespace models {

    v_model::v_model()
    : mem_dim(dim2(16,16)), work_dim(dim2(16,16)), item_dim(dim2(16,16)) // to revert to 128,128
    {
    }

    v_model::~v_model(){
    }

    void v_model::add_revision(imodel::object* obj){
        static __a_timer time("add_revision");
        static __a_timer time2("new_layout");
        static __a_timer time3("map insert");
        static __a_timer time4(">> operator for layout");
        time2.begin();
        v_model::layout* l = new v_model::layout(obj->get_dim(), obj->get_t_size());
        time2.end();
        time3.begin();
        *const_cast<size_t*>(&l->sid) = this->map.insert(l);
        *const_cast<const size_t **>(&l->gid) = channel.id().first;
        time3.end();
        time4.begin();
        (*l) >> this->mem_dim, this->work_dim, this->item_dim;
        time4.end();
        time.begin();
        obj->add_revision(l);
        time.end();
    }

    void v_model::update_revision(imodel::revision* r, channels::group* placement){

    }

    v_model::revision* v_model::get_revision(size_t* hash, size_t hash_len, size_t id) const {
        return ((layout*)this->map.find(id))->revision;
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

    void* solidify(const imodel::object& o){
        imodel::revision& r = o.revision(0);
        imodel::layout& l = r.get_layout();
        size_t iterator = 0;
        char* memory = NULL;

        for(size_t j=0; j < l.get_mem_grid_dim().x; j++)
            for(size_t jj=0; jj < l.get_mem_dim().x; jj++)
                for(size_t i=0; i < l.get_mem_grid_dim().y; i++){
                    if(r.block(i,j)->valid()){
                        memory = (char*)realloc(memory, (iterator+1)*l.get_mem_lda());
                        memcpy(memory+iterator*l.get_mem_lda(),                         // copy to
                               &((char*)r(i,j))[jj*l.get_mem_lda()],                    // copy from
                               l.get_mem_lda());                                        // of size
                        iterator++;
                    }
                }
        return memory;
    }

    void disperse(void* data, imodel::object& o){
        imodel::revision& current = o.revision(0);
        imodel::revision& updated = o.revision(1);
        imodel::layout& l = current.get_layout();
        size_t iterator = 0;
        char* memory = (char*)data;

        for(size_t j=0; j < l.get_mem_grid_dim().x; j++)
            for(size_t jj=0; jj < l.get_mem_dim().x; jj++)
                for(size_t i=0; i < l.get_mem_grid_dim().y; i++){
                    if(current.block(i,j)->valid()){
                        memcpy(&((char*)updated(i,j))[jj*l.get_mem_lda()],  // copy from
                               memory+iterator*l.get_mem_lda(),             // copy to
                               l.get_mem_lda());                            // of size
                        iterator++;
                    }
                }
        free(data);
    }

    // }}}

} }
