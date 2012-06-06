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

    inline revision& model::add_revision(history* o){
        dim2 block_dim;
        if(this->mem_dim > o->get_cached_dim().max()) block_dim = o->get_cached_dim();
        else block_dim = this->mem_dim;
        layout* l = new layout(o->get_t_size(), block_dim, o->get_cached_dim());
        this->insert(l);
        return o->add_state(l);
    }

    inline bool model::is_atomic(const history* o){
        if(o->back() == NULL) return (this->mem_dim > o->get_cached_dim().max()); 
        return (o->back()->content->grid_dim == 1);
    }

    inline size_t model::get_block_lda(history* o){
        return o->back()->content->mem_dim.y;
    }

    inline dim2 model::get_current_dim(const history* o){
        return o->get_cached_dim();
    }

    inline void model::set_current_dim(history* o, dim2 dim){
        if(o->back() != NULL)
            o->back()->set_dim(dim);
        o->cache_dim(dim);
    }

    inline size_t model::time(const history* o){
        if(o->back() == NULL) 
            this->add_revision(const_cast<history*>(o));
        return o->time();
    }

    inline layout* model::get_layout(size_t id) const {
        return (layout*)this->map.find(id);
    }

    inline model& model::operator>>(dim2 mem_dim){
        this->mem_dim  = mem_dim;
        return *this;
    }

} } }
