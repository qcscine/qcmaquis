#include "ambient/ambient.h"
#include "ambient/models/v_model.h"

#include "ambient/models/operation/operation.h"
    
namespace ambient { namespace models {

    v_model::revision::revision(imodel::layout* l)
    : layout(l), initialization(NULL), reduction(NULL)
    {
        this->layout->set_revision(this);
    };

    v_model::revision::~revision(){
        delete this->layout;
    }

    void v_model::revision::set_dim(dim2 dim){
        this->layout->set_dim(dim);
    }

    std::pair<size_t*,size_t> v_model::revision::id(){
        return this->layout->id();
    }

    imodel::layout& v_model::revision::get_layout(){
        return *this->layout;
    }

    imodel::layout::entry* v_model::revision::block(size_t i, size_t j){
        return this->layout->get(i, j);
    }

    imodel::layout::entry& v_model::revision::operator()(size_t i, size_t j){
        return *(imodel::layout::entry*)this->block(i,j);
    }

    void v_model::revision::add_modifier(models::imodel::modifier* m){
        this->modifiers.push_back(m);
    }

    channels::group* v_model::revision::get_placement(){
        return this->placement;
    }

    void v_model::revision::reduce(void(*fp)(void*,void*)){
       this->reduction = fp;
    }

    void v_model::revision::init(void(*fp)()){
       this->initialization = fp;
    }

} }
