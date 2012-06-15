namespace ambient { namespace controllers { namespace velvet {

    inline layout& revision_sub::get_layout()  { return ((revision*)this)->get_layout(); }
    inline revision* revision_sub::get_parent(){ return ((revision*)this)->get_parent(); }

    inline layout::entry& c_revision::operator()(size_t x, size_t y){
        layout::entry& e = ((revision*)this)->block(x,y);
        if(!e.valid()) ambient::controller.alloc_block(this->get_layout().spec, this->get_layout(), x, y);
        return e;
    }
    inline layout::entry& p_revision::operator()(size_t x, size_t y){
        layout::entry& e = ((revision*)this)->block(x,y);
        if(!e.valid()) ambient::controller.calloc_block(this->get_layout().spec, this->get_layout(), x, y);
        return e;
    }
    inline layout::entry& w_revision::operator()(size_t x, size_t y){
        layout::entry& e = ((revision*)this)->block(x,y);
        if(!e.valid()){
            layout::entry& parent = this->get_parent()->block(x,y);
            if(parent.occupied()) ambient::controller.alloc_block(this->get_layout().spec, this->get_layout(), x, y);
            else e.swap(parent);
        }
        return e;
    }
    inline layout::entry& r_revision::operator()(size_t x, size_t y){
        layout::entry& e = ((revision*)this)->block(x,y);
        if(!e.valid()){
            layout::entry& parent = this->get_parent()->block(x,y);
            if(parent.occupied()) ambient::controller.alloc_block(this->get_parent()->get_layout().spec, this->get_layout(), x, y);
            else e.swap(parent);
        }
        return e;
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
    inline size_t iteratable<T>::get_thread_revision_base() const {
        return this->thread_revision_base[ctxt.get_tid()];
    }

    template<class T>
    inline void iteratable<T>::set_thread_revision_base(size_t r){
        this->thread_revision_base[ctxt.get_tid()] = r;
    }

} } }
