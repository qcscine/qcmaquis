#include "ambient/ambient.h"
#include "ambient/models/v_model.h"

#include "ambient/models/operation/operation.h"
    
namespace ambient { namespace models {

    v_model::object::object()
    : placement(NULL) 
    {
        this->base = 0;
    }

    v_model::object::~object(){
        for(size_t i=0; i < this->revisions.size(); i++)
            delete this->revisions[i];
    }

    v_model::revision& v_model::object::revision(size_t offset) const {
        if(this->revisions.size() == 0) ambient::model.add_revision(const_cast<v_model::object*>(this));
        return *this->revisions[this->base + offset];
    }

    void v_model::object::add_revision(imodel::layout* l){
        this->revisions.push_back(new v_model::revision(l));
    }

    dim2 v_model::object::get_dim() const {
        return this->dim;
    }

    size_t v_model::object::get_t_size() const {
        return this->t_size;
    }

    channels::group* v_model::object::get_placement(){
        return this->placement;
    }

} }
