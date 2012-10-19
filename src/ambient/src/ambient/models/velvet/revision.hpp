namespace ambient { namespace models { namespace velvet {

    inline void* revision::operator new (size_t size){
        return ambient::static_memory::malloc<revision>(); 
    }

    inline void revision::operator delete (void* ptr){
        ambient::static_memory::free<revision>(ptr);
    }

    inline revision::revision(memspec* spec, void* g)
    : spec(spec), clean(g == NULL), generator(g), header(NULL), data(NULL), users(0) 
    {
    }

    inline void revision::embed(void* memory, size_t bound){
        this->header = memory;
        this->data = (void*)((size_t)memory + bound);
    }

    inline void revision::reset_generator(){
        this->generator = NULL;
    }

    inline void* revision::get_generator(){
        return this->generator;
    }

    inline void revision::swap(revision& e){
        this->data = e.data;
        this->header = e.header;
        e.header = NULL;
    }

    inline void* revision::get_memory(){
        return this->header;
    }

    inline bool revision::valid(){
        return (this->data != NULL);
    }

    inline bool revision::occupied(){
        return (this->users != 0);
    }

    inline void revision::release(){
        --this->users;
    }

    inline void revision::use(){
        ++this->users;
    }

    // }}}

} } }
