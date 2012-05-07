#include "ambient/ambient.h"
#include "ambient/models/v_model.h"

namespace ambient { namespace models {

    v_model::object::object()
    : revision_base(0)
    {
        this->thread_revision_base = (size_t*)calloc(controller.get_num_threads(), sizeof(size_t));
    }

    v_model::object::~object(){
        for(size_t i=0; i < this->revisions.size(); i++)
            delete this->revisions[i];
    }

    v_model::revision& v_model::object::revision(size_t offset) const {
        if(this->revisions.size() == 0) ambient::model.add_revision(const_cast<v_model::object*>(this));
        assert(ctxt.get_revision_base(this) + offset < this->revisions.size());

        return *this->revisions[ctxt.get_revision_base(this) + offset];
    }

    void v_model::object::add_revision(imodel::layout* l){
        if(this->revisions.size() == 0) l->mark(0, 0);
        this->set_revision_base(this->revisions.size());
        this->revisions.push_back(new v_model::revision(this, l));
    }

    dim2 v_model::object::get_dim() const {
        if(this->revisions.size() == 0) return this->init_dim;
        return this->revision(0).get_dim();
    }

    void v_model::object::set_dim(dim2 dim){
        if(this->revisions.size() == 0) this->init_dim = dim;
        else this->revision(0).set_dim(dim); // used in pt_set_dim after pushes
    }

    size_t v_model::object::get_t_size() const {
        return this->t_size;
    }

    size_t v_model::object::get_revision_base() const {
        return this->revision_base;
    }

    size_t v_model::object::get_thread_revision_base() const {
        return this->thread_revision_base[ctxt.get_tid()];
    }

    void v_model::object::set_revision_base(size_t r){
        this->revision_base = r;
    }

    void v_model::object::set_thread_revision_base(size_t r){
        this->thread_revision_base[ctxt.get_tid()] = r;
    }

} }
