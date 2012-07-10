namespace ambient { namespace models { namespace velvet {

    inline model::model()
    : mem_dim(dim2(10240,10240))
    {
    }

    inline model::~model(){
    }

    inline void model::insert(revision* r){
        *const_cast<size_t*>(&r->sid) = this->map.insert(r);
    }

    template<typename T>
    inline revision* model::init_revision(T* o){
        //if(this->mem_dim < o->spec.dim.max()) block = this->mem_dim; else 
        dim2 block = o->spec.dim;
        o->spec.latch<T::value_type>(block);
        revision* r = new revision(&o->spec, true);
        o->add_state(r);
        return r;
    }

    inline revision* model::add_revision(history* o){
        revision* r = new revision(&o->spec);
        o->add_state(r);
        return r;
    }

    inline bool model::is_atomic(const history* o){
        if(o->back() == NULL) return (this->mem_dim > o->spec.dim.max()); 
        return (o->spec.grid == 1);
    }

    template<typename T>
    inline size_t model::time(const T* o){
        if(o->back() == NULL)
            this->init_revision(const_cast<T*>(o));
        return o->time();
    }

    inline revision* model::get_revision(size_t id) const {
        return (revision*)this->map.find(id);
    }

    inline model& model::operator>>(dim2 mem_dim){
        this->mem_dim  = mem_dim;
        return *this;
    }

} } }
