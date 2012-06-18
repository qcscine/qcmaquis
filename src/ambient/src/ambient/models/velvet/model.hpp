namespace ambient { namespace models { namespace velvet {

    inline model::model()
    : mem_dim(dim2(1024,1024))
    {
    }

    inline model::~model(){
    }

    inline void model::insert(layout* l){
        *const_cast<size_t*>(&l->sid) = this->map.insert(l);
    }

    template<typename T>
    inline revision& model::init_revision(T* o){
        dim2 block;
        if(this->mem_dim < o->get_cached_dim().max()) block = this->mem_dim;
        else block = o->get_cached_dim();
        o->init_spec = new memspec(block.square()*sizeof(typename T::value_type), 
                                   block, o->get_cached_dim());
        o->spec = new memspec(*o->init_spec);
        layout* l = new layout(o->init_spec);
        o->init_spec->clean = true;
        //this->insert(l);
        return o->add_state(l);
    }

    inline revision& model::add_revision(history* o){
        layout* l = new layout(o->spec);
        //this->insert(l);
        return o->add_state(l);
    }

    inline bool model::is_atomic(const history* o){
        if(o->back() == NULL) return (this->mem_dim > o->get_cached_dim().max()); 
        return (o->spec->grid == 1);
    }

    inline size_t model::get_block_lda(history* o){
        return o->spec->grid.y;
    }

    inline dim2 model::get_current_dim(const history* o){
        return o->get_cached_dim();
    }

    inline void model::set_current_dim(history* o, dim2 dim){
        o->cache_dim(dim);
    }

    template<typename T>
    inline size_t model::time(const T* o){
        if(o->back() == NULL)
            this->init_revision(const_cast<T*>(o));
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
