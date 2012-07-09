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

    inline layout::layout(memspec* spec, bool clean)
    : spec(spec), clean(clean), lda(spec->grid.y), entries(spec->grid.square())
    {
        std::for_each(entries.begin(), entries.end(), constructor());
    }

    inline void layout::embed(void* memory, size_t x, size_t y, size_t bound){
        this->get(x,y).set_memory(memory, bound);
    }

    inline layout::entry& layout::get(size_t x, size_t y){
        return *this->entries[x*lda+y];
    }

    //inline size_t layout::get_master(){
    //    return this->master;
    //}

    // }}}

    // {{{ layout::entry //

    //inline bool layout::entry::trylock(){
    //    if(this->locked) return false;
    //    this->locked = true;
    //    return true;
    //}

    //inline void layout::entry::unlock(){
    //    this->locked = false;
    //}

    inline layout::entry::entry()
    : header(NULL), data(NULL)//, request(false), locked(false)
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

    //inline bool layout::entry::requested(){
    //    return this->request;
    //}

    inline std::list<controllers::velvet::cfunctor*>& layout::entry::get_assignments(){
        return this->assignments;
    }

    //inline std::list<size_t>& layout::entry::get_path(){
    //    return this->path;
    //}

    // }}}

} } }
