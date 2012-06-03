namespace ambient { namespace models { namespace velvet {

    inline model::model()
    : mem_dim(dim2(256,256))
    {
    }

    inline model::~model(){
    }

    inline void model::insert(layout* l){
        *const_cast<size_t*>(&l->sid) = this->map.insert(l);
    }

    inline void model::add_revision(history* obj){
        layout* l = new layout(obj->get_t_size());
        l->set_dim(obj->get_dim());
        this->insert(l);
        obj->add_state(l);
    }

    inline void model::update_revision(revision* r, group* placement){

    }

    inline layout* model::get_layout(size_t id) const {
        return (layout*)this->map.find(id);
    }

    inline model& model::operator>>(dim2 mem_dim){
        this->mem_dim  = mem_dim;
        return *this;
    }

} } }
