namespace ambient { namespace models { namespace velvet {

    inline revision::revision(memspec* spec, bool clean)
    : spec(spec), clean(clean), generator(NULL) { 
    }

    inline revision::~revision(){
        free(content.header);
    }

    inline void revision::embed(void* memory, size_t x, size_t y, size_t bound){
        this->content.set_memory(memory, bound);
    }

    inline revision::entry& revision::block(size_t x, size_t y){
        return this->content;
    }

    inline void revision::add_modifier(sfunctor* m){
        this->modifiers.push_back(m);
    }

    inline std::list<sfunctor*>& revision::get_modifiers(){
        return this->modifiers;
    }

    inline void revision::set_generator(sfunctor* m){
        this->generator = m;
    }

    inline void revision::reset_generator(){
        this->generator = NULL;
    }

    inline sfunctor* revision::get_generator(){
        return this->generator;
    }

    // {{{ revision::entry //

    inline bool revision::entry::trylock(){
        //if(this->locked) return false;
        //this->locked = true;
        return true;
    }

    inline void revision::entry::unlock(){
        //this->locked = false;
    }

    inline revision::entry::entry()
    : header(NULL), data(NULL){ }

    inline void revision::entry::swap(entry& e){
        this->data = e.data;
        this->header = e.header;
        e.header = NULL;
    }

    inline void revision::entry::set_memory(void* memory, size_t bound){
        this->header = memory;
        this->data = (void*)((size_t)memory + bound);
    }

    inline void* revision::entry::get_memory(){
        return this->header;
    }

    inline bool revision::entry::valid(){
        return (this->data != NULL);
    }

    inline bool revision::entry::occupied(){
        return (this->data == NULL);
    }

    inline std::list<controllers::velvet::cfunctor*>& revision::entry::get_assignments(){
        return this->assignments;
    }

    // }}}

} } }
