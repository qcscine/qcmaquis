namespace ambient { namespace controllers { namespace velvet {

    inline revision* revision_sub::get_parent(){ return ((revision*)this)->get_parent(); }

    inline revision::entry& c_revision::operator()(size_t x, size_t y){
        revision::entry& e = ((revision*)this)->content;
        if(!e.valid()){
            if(((revision*)this)->clean) ambient::controller.calloc_block(*(revision*)this, x, y);
            else ambient::controller.alloc_block(*(revision*)this, x, y);
        }
        return e;
    }
    inline revision::entry& p_revision::operator()(size_t x, size_t y){
        revision::entry& e = ((revision*)this)->content;
        if(!e.valid()) ambient::controller.calloc_block(*(revision*)this, x, y);
        return e;
    }
    inline revision::entry& w_revision::operator()(size_t x, size_t y){
        revision::entry& e = ((revision*)this)->content;
        if(!e.valid()){
            revision::entry& parent = this->get_parent()->content;
            if(parent.occupied()) ambient::controller.alloc_block(*(revision*)this, x, y);
            else e.swap(parent);
        }
        return e;
    }
    inline revision::entry& r_revision::operator()(size_t x, size_t y){
        revision::entry& e = ((revision*)this)->content;
        if(!e.valid()){
            revision::entry& parent = this->get_parent()->content;
            if(parent.occupied()){
                if(this->get_parent()->clean) ambient::controller.calloc_block(*(revision*)this, x, y);
                else ambient::controller.alloc_block(*(revision*)this, x, y);
            }else e.swap(parent);
        }
        return e;
    }
    
    template<class T>
    inline iteratable<T>::iteratable(dim2 dim) : T(dim) { }

    template<class T>
    inline size_t iteratable<T>::get_thread_revision_base() const {
        return this->thread_revision_base[ctxt.get_tid()];
    }

    template<class T>
    inline void iteratable<T>::set_thread_revision_base(size_t r){
        this->thread_revision_base[ctxt.get_tid()] = r;
    }

} } }
