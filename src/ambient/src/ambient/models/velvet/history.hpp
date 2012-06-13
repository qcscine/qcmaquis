#define NUM_THREADS 1 // controller.get_num_threads()

namespace ambient { namespace models { namespace velvet {

    inline history::history()
    : current(NULL), start(0)
    {
    }

    inline history::~history(){
        size_t size = this->content.size();
        for(size_t i = start; i < size; i++)
            delete this->content[i];
    }

    inline void history::clean(){
        while(this->content[start] != this->current)
            delete this->content[start++];
    }

    inline revision& history::add_state(layout* l){
        revision* parent = this->current;
        this->current = new revision(l);
        this->current->parent = parent;
        this->content.push_back(this->current);
        return *this->current;
    }

    inline revision* history::back() const {
        return this->current;
    }

    inline size_t history::time() const {
        return this->content.size()-1;
    }

    inline dim2 history::get_cached_dim() const {
        return this->dim;
    }

    inline void history::cache_dim(dim2 dim){
        this->dim = dim; 
    }

} } }
