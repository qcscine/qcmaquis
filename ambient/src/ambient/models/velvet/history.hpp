namespace ambient { namespace models { namespace velvet {

    inline void* history::operator new (size_t size){
        return ambient::static_memory::malloc<history>();
    }

    inline void history::operator delete (void* ptr){
        ambient::static_memory::free<history>(ptr);
    }

    inline history::history(dim2 dim, size_t ts) : current(NULL), spec(dim, ts) { 
        this->content.reserve(2); 
    }

    inline revision* history::add_state(void* g){
        revision* r = new revision(&spec, g); 
        r->parent = this->current;
        this->content.push_back(r);
        this->current = r;
        return r;
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
