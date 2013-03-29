namespace ambient { namespace models { namespace velvet {

    inline void model::index(revision* r){
        r->sid = this->sid++;
        this->sid %= MAX_SID;
    }

    inline void model::index(transformable* v){
        v->sid = this->sid++;
        this->sid %= MAX_SID;
    }

    template<ambient::locality L>
    inline void model::add_revision(history* o, void* g){
        o->add_state<L>(g);
    }

    inline void model::use_revision(history* o){
        o->back()->use();
    }

    inline void model::touch(const history* o){
        if(o->back() == NULL) 
            this->add_revision<ambient::common>(const_cast<history*>(o));
    }

    inline size_t model::time(const history* o){
        this->touch(o);
        return o->time();
    }

    inline bool model::feeds(const revision* r){
        return (r->state == ambient::local);
    }

    inline bool model::common(const revision* r){
        return (r->state == ambient::common);
    }

} } }
