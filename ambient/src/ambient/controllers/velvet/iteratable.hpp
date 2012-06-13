namespace ambient { namespace controllers { namespace velvet {


    inline layout::entry& check_revision::operator()(size_t x, size_t y){
        layout::entry& e = ((revision*)this)->block(x,y);
        if(!e.valid()) return ambient::controller.alloc_block(this->get_layout(), x, y);
        return e;
    }

    inline layout& check_revision::get_layout(){
        return ((revision*)this)->get_layout();
    }

    inline layout::entry& reuse_revision::operator()(size_t x, size_t y){
        layout::entry& e = ((revision*)this)->block(x,y);
        if(!e.valid()){
            layout::entry& parent = this->get_parent()->block(x,y);
            if(parent.occupied()) ambient::controller.alloc_block(this->get_layout(), x, y);
            else e.swap(parent);
        }
        return e;
    }

    inline layout& reuse_revision::get_layout(){
        return ((revision*)this)->get_layout();
    }

    inline revision* reuse_revision::get_parent(){
        return ((revision*)this)->get_parent();
    }
    
    template<class T>
    inline iteratable<T>::iteratable(){
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
