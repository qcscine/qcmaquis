#include "ambient/ambient.h"
#include "ambient/models/v_model.h"

namespace ambient { namespace models {

    void null_i(models::v_model::object& a){
        dim2 idx = ctxt.get_block_id();
        dim2 dim = a.revision(0).get_layout().get_mem_dim();
        memset((double*)a.revision(0)(idx.y,idx.x), 0, dim.y*dim.x*a.get_t_size());
    }

    v_model::revision::revision(imodel::object* o, imodel::layout* l)
    : object(o), layout(l), initialization((voidfp)null_i), reduction(NULL), placement(NULL), generator(NULL)
    {
        this->number = o->get_revision_base(); // debug
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

    imodel::object& v_model::revision::get_object(){
        return *this->object;
    }

    imodel::layout& v_model::revision::get_layout(){
        return *this->layout;
    }

    imodel::layout::entry* v_model::revision::block(size_t i, size_t j){
        return this->layout->get(i, j);
    }

    imodel::layout::entry& v_model::revision::operator()(size_t i, size_t j){
        return controller.ufetch_block(*this, i, j);
    }

    void v_model::revision::add_modifier(models::imodel::modifier* m){
        this->modifiers.push_back(m);
    }

    std::list<imodel::modifier*>& v_model::revision::get_modifiers(){
        return this->modifiers;
    }

    channels::group* v_model::revision::get_placement(){
        return this->placement;
    }

    void v_model::revision::set_placement(channels::group* grp){
        this->placement = grp;
    }

    void v_model::revision::set_generator(imodel::modifier* m){
        this->generator = m;
    }

    imodel::modifier* v_model::revision::get_generator(){
        return this->generator;
    }

    void v_model::revision::reduce(void(*fp)()){
       this->reduction = fp;
    }

    void v_model::revision::init(void(*fp)()){
       this->initialization = fp;
    }

    v_model::revision::voidfp v_model::revision::get_reduce(){
        return this->reduction;
    }

     v_model::revision::voidfp v_model::revision::get_init(){
        return this->initialization;
    }

} }
