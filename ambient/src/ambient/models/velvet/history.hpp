#define NUM_THREADS 1 // controller.get_num_threads()

namespace ambient { namespace models { namespace velvet {

    inline history::history(size_t ts)
    : t_size(ts)
    {
        layout* l = new layout(ts);
        ambient::model.insert(l);
        this->add_state(l);
    }

    inline history::~history(){
        for(size_t i=0; i < this->content.size(); i++)
            delete this->content[i];
    }

    inline void history::add_state(layout* l){
        this->current = new revision(l);
        this->content.push_back(this->current);
    }

    inline revision* history::back(){
        return this->current;
    }

    inline dim2 history::get_dim() const {
        return this->current->get_dim();
    }

    inline void history::set_dim(dim2 dim){
        this->current->set_dim(dim); // used in pt_set_dim after pushes
    }

    inline size_t history::get_t_size() const {
        return this->t_size;
    }

    inline size_t history::time() const {
        return this->content.size()-1;
    }

} } }
