namespace ambient { namespace models { namespace velvet {

    inline void model::insert(revision* r){
        *const_cast<size_t*>(&r->sid) = this->map.insert(r);
    }

    inline revision* model::add_revision(history* o, void* g){
        return o->add_state(g);
    }

    inline void model::use_revision(history* o){
        o->back()->use();
    }

    inline size_t model::time(const history* o){
        if(o->back() == NULL)
            this->add_revision(const_cast<history*>(o));
        return o->time();
    }

    inline revision* model::get_revision(size_t id) const {
        return (revision*)this->map.find(id);
    }

} } }
