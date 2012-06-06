#define NUM_THREADS 1 // controller.get_num_threads()

namespace ambient { namespace models { namespace velvet {

    inline history::history(size_t ts)
    : t_size(ts), current(NULL)
    {
    }

    inline history::~history(){
        for(size_t i=0; i < this->content.size(); i++)
            delete this->content[i];
    }

    inline revision& history::add_state(layout* l){
        this->current = new revision(l);
        this->content.push_back(this->current);
        return *this->current;
    }

    inline revision* history::back() const {
        return this->current;
    }

    inline size_t history::get_t_size() const {
        return this->t_size;
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
