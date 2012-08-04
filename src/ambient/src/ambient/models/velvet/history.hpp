namespace ambient { namespace models { namespace velvet {

    inline void* history::operator new (size_t size){
        return boost::singleton_pool<ambient::utils::empty, sizeof(history)>::malloc(); 
    }

    inline void history::operator delete (void* ptr){
        boost::singleton_pool<ambient::utils::empty, sizeof(history)>::free(ptr); 
    }

    inline history::history(dim2 dim, size_t ts) : current(NULL), start(0), spec(dim, ts), references(0) { 
        this->content.reserve(2); 
    }

    inline history::~history(){
        size_t size = this->content.size();
        for(size_t i = start; i < size; i++)
            delete this->content[i];
    }

    inline void history::clean(){
        //if(current != NULL)
        //while(this->content[start] != this->current){
        //    if(this->content[start]->content.get_assignments().size() != 0) break;
        //    else printf("Going to delete the revision! (%lu)\n", this->content[start]->content.get_assignments().size());
        //    delete this->content[start++];
        //}
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
