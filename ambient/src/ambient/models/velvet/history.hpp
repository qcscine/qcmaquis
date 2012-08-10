namespace ambient { namespace models { namespace velvet {

    inline void* history::operator new (size_t size){
        return boost::singleton_pool<ambient::utils::empty, sizeof(history)>::malloc(); 
    }

    inline void history::operator delete (void* ptr){
        boost::singleton_pool<ambient::utils::empty, sizeof(history)>::free(ptr); 
    }

    inline history::history(dim2 dim, size_t ts) : current(NULL), spec(dim, ts) { 
        this->content.reserve(2); 
    }

    inline history::~history(){
        size_t size = this->content.size();
        for(size_t i = 0; i < size; i++){
            free(this->content[i]->header);
            delete this->content[i];
        }
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
