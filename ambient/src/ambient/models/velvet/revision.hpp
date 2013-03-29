namespace ambient { namespace models { namespace velvet {

    inline void* revision::operator new (size_t size){
        return ambient::pool.malloc<revision>(); 
    }

    inline void revision::operator delete (void* ptr){
        ambient::pool.free<revision>(ptr);
    }

    inline revision::revision(size_t extent, void* g, ambient::locality l)
    : extent(extent), generator(g), state(l), 
      data(NULL), transfer(NULL), users(0), 
      region(DEFAULT_REGION)
    {
    }

    inline void revision::embed(void* memory){
        data = memory;
    }

    inline void revision::reuse(revision& r){
        data     = r.data;
        region   = r.region;
        r.region = DELEGATED_REGION;
    }

    inline void revision::use(){
        ++users;
    }

    inline void revision::release(){
        --users;
    }

    inline void revision::complete(){
        generator = NULL;
    }

    inline bool revision::locked(){
        return (users != 0);
    }

    inline bool revision::valid(){
        return (data != NULL);
    }

    // }}}

} } }
