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

    inline history::~history(){
        int size = this->content.size();
        for(int i = 0; i < size; i++) spec.free(this->content[i]->header);
        for(int i = 0; i < size; i++) ambient::static_memory::free<revision>(this->content[i]); 
    }

    inline void history::add_state(revision* r){
        r->parent = this->current;
        this->content.push_back(r);
        this->current = r;
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
