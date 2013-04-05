namespace ambient { namespace models { namespace velvet {

    inline void* history::operator new (size_t size){
        return ambient::pool.malloc<history>();
    }

    inline void history::operator delete (void* ptr){
        ambient::pool.free<history>(ptr);
    }

    inline history::history(dim2 dim, size_t ts) : current(NULL), spec(dim, ts) {
        this->clock = ambient::model.clock;
        this->content.reserve(2); 
    }

    template<ambient::locality L>
    inline void history::add_state(void* g){
        revision* r = new revision(spec.size, g, L); 
        this->content.push_back(r);
        this->current = r;
        ambient::model.index(r);
    }

    inline void history::fuse(const history* src){
        if(src->weak()) return;
        revision* r = src->back();
        this->content.push_back(r);
        this->current = r;
        // do not deallocate or reuse
        if(!r->valid()) r->region++;
        assert(!r->valid() || r->region > 0); // can't rely on bulk memory
        r->use();
    }
        
    inline size_t history::time() const {
        return this->content.size()-1;
    } 

    inline revision* history::back() const {
        return this->current;
    }

    inline bool history::weak() const {
        return (this->back() == NULL);
    }

} } }
