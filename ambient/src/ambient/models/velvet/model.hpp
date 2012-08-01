namespace ambient { namespace models { namespace velvet {

    inline void model::insert(revision* r){
        *const_cast<size_t*>(&r->sid) = this->map.insert(r);
    }

    template<typename T>
    inline revision* model::init_revision(T* o){
        revision* r = new revision(&o->spec, true);
        o->add_state(r);
        return r;
    }

    inline revision* model::add_revision(history* o){
        revision* r = new revision(&o->spec);
        o->add_state(r);
        return r;
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

} } }
