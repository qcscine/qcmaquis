#define NUM_THREADS 1 // controller.get_num_threads()

namespace ambient { namespace models { namespace velvet {

    inline history::history(dim2 dim) : current(NULL), start(0), spec(dim) { }

    inline history::~history(){
        size_t size = this->content.size();
        for(size_t i = start; i < size; i++)
            delete this->content[i];
    }

    inline void history::clean(){
        while(this->content[start] != this->current)
            delete this->content[start++];
    }

    inline revision& history::add_state(revision* r){
        r->parent = this->current;
        this->content.push_back(r);
        this->current = r;
        return *r;
    }

    inline revision* history::back() const {
        return this->current;
    }

    inline size_t history::time() const {
        return this->content.size()-1;
    }

    inline bool history::weak() const {
        return (this->back() == NULL);
    }

} } }
