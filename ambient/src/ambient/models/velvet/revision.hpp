namespace ambient { namespace models { namespace velvet {

    inline void* revision::operator new (size_t size){
        return ambient::pool.malloc<revision>(); 
    }

    inline void revision::operator delete (void* ptr){
        ambient::pool.free<revision>(ptr);
    }

    inline revision::revision(size_t extent, void* g)
    : extent(extent), generator(g), users(0)
    {
        state = (generator == NULL) ? PURE : VOID;
        header = data = NULL;
    }

    inline void revision::embed(void* memory, size_t bound){
        this->header = memory;
        this->data = (void*)((size_t)memory + bound);
    }

    inline void revision::reuse(revision& r){
        this->data   = r.data;
        this->header = r.header;
        r.header     = NULL;
    }

    inline void revision::use(){
        ++this->users;
    }

    inline void revision::release(){
        --this->users;
    }

    inline void revision::complete(){
        this->generator = NULL;
    }

    inline bool revision::locked(){
        return (this->users != 0);
    }

    inline bool revision::remote(){
        return (this->state == VOID);
    }

    inline bool revision::origin(){
        return (this->state == ORIG);
    }

    inline bool revision::valid(){
        return (this->data != NULL); //(this->state == ORIG || 
                                     // this->state == COPY );
    }

    // }}}

} } }
