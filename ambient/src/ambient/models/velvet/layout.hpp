#include "ambient/utils/ceil.h"
//#define LAYOUT_ACCESS_CHECK


namespace ambient { namespace models { namespace velvet {

    // {{{ layout model

    struct deleter {
       inline void operator()(layout::entry* e){ free(e->header); delete e; }
    };
    struct constructor {
       inline void operator()(layout::entry*& e){ e = new layout::entry(); }
    };

    inline layout::~layout(){
        std::for_each(entries.begin(), entries.end(), deleter());
    }

    inline layout::layout(size_t t_size, dim2 b_size, dim2 size)
    : spec(t_size*b_size.square()), grid_dim(__a_ceil(size.x / b_size.x), __a_ceil(size.y / b_size.y)), 
      dim(size), mem_dim(b_size), entries(grid_dim.square()) 
    {
        std::for_each(entries.begin(), entries.end(), constructor());
    }

    inline const memspec& layout::get_spec() const {
        return this->spec;
    }

    inline void layout::embed(void* memory, size_t x, size_t y, size_t bound){
        this->get(x,y).set_memory(memory, bound);
    }

    inline layout::entry& layout::get(size_t x, size_t y){
#ifdef LAYOUT_ACCESS_CHECK
        if(y >= this->get_grid_dim().y || x >= this->get_grid_dim().x)
        printf("%ld: Trying to access %ld x %ld of %ld x %ld\n", this->sid, x, y, this->get_grid_dim().x, this->get_grid_dim().y);
#endif
        return *this->entries[x*grid_dim.y+y];
    }

    inline size_t layout::get_master(){
        return this->master;
    }

    inline dim2 layout::get_dim() const {
        return this->dim;
    }

    inline dim2 layout::get_mem_dim() const { 
        return this->mem_dim;
    }

    inline dim2 layout::get_grid_dim() const {
        return this->grid_dim;
    }

    // }}}

    // {{{ layout::entry //

    inline bool layout::entry::trylock(){
        if(this->locked) return false;
        this->locked = true;
        return true;
    }

    inline void layout::entry::unlock(){
        this->locked = false;
    }

    inline layout::entry::entry()
    : header(NULL), data(NULL), request(false), locked(false)
    {
    }

    inline void layout::entry::swap(entry& e){
        this->data = e.data;
        this->header = e.header;
        e.header = NULL;
    }

    inline void layout::entry::set_memory(void* memory, size_t bound){
        this->header = memory;
        this->data = (void*)((size_t)memory + bound);
    }

    inline void* layout::entry::get_memory(){
        return this->header;
    }

    inline bool layout::entry::valid(){
        return (this->data != NULL);
    }

    inline bool layout::entry::occupied(){
        return (this->data == NULL);
    }

    inline bool layout::entry::requested(){
        return this->request;
    }

    inline std::list<controllers::velvet::cfunctor*>& layout::entry::get_assignments(){
        return this->assignments;
    }

    inline std::list<size_t>& layout::entry::get_path(){
        return this->path;
    }

    // }}}

} } }
