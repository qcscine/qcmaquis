namespace ambient { namespace controllers { namespace velvet {


    inline layout::entry& slow_revision::operator()(size_t x, size_t y){
        layout::entry* e = ((revision*)this)->block(x,y);
        if(e->header == NULL) return ambient::controller.alloc_block(this->get_layout(), x, y);
        return *e;
    }

    inline layout& slow_revision::get_layout(){
        return ((revision*)this)->get_layout();
    }
    
    inline layout::entry& fast_revision::operator()(size_t x, size_t y){
        // while(!this->valid()) pthread_yield(); // can be done if MPI enabled
        return *((revision*)this)->block(x,y);
    }

    inline layout& fast_revision::get_layout(){
        return ((revision*)this)->get_layout();
    }

    template<class T>
    inline iteratable<T>::iteratable(size_t ts)
    : T(ts) {
        this->thread_revision_base = (size_t*)calloc(NUM_THREADS, sizeof(size_t));
    }

    template<class T>
    inline iteratable<T>::~iteratable(){
        free(this->thread_revision_base);
    }

    template<class T>
    inline revision& iteratable<T>::state(size_t offset) const {
        //assert(ctxt.get_revision_base(this) + offset < ((T*)this)->content.size());
        return *((T*)this)->content[ctxt.get_revision_base(this) + offset];
    } 

    template<class T>
    inline size_t iteratable<T>::get_thread_revision_base() const {
        return this->thread_revision_base[ctxt.get_tid()];
    }

    template<class T>
    inline void iteratable<T>::set_thread_revision_base(size_t r){
        this->thread_revision_base[ctxt.get_tid()] = r;
    }

} } }
