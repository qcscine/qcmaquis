namespace ambient { namespace models { namespace velvet {

    inline void* revision::operator new (size_t size){
        return boost::singleton_pool<ambient::utils::empty, sizeof(revision)>::malloc(); 
    }

    inline void revision::operator delete (void* ptr){
        boost::singleton_pool<ambient::utils::empty, sizeof(revision)>::free(ptr); 
    }

    inline revision::revision(memspec* spec, bool clean)
    : spec(spec), clean(clean), generator(NULL), header(NULL), data(NULL){
    }

    inline void revision::embed(void* memory, size_t bound){
        this->header = memory;
        this->data = (void*)((size_t)memory + bound);
    }

    inline void revision::set_generator(void* m){
        this->generator = m;
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
        return (this->data == NULL);
    }

    // }}}

} } }
