namespace ambient { namespace controllers { namespace velvet {

    inline revision* revision_sub::get_parent(){ return ((revision*)this)->get_parent(); }

    template<typename T>
    inline c_revision::operator T* (){
        revision::entry& e = ((revision*)this)->content;
        if(!e.valid()){
            if(((revision*)this)->clean) ambient::controller.calloc(*(revision*)this);
            else ambient::controller.alloc(*(revision*)this);
        }
        return (T*)e;
    }
    template<typename T>
    inline p_revision::operator T* (){
        revision::entry& e = ((revision*)this)->content;
        if(!e.valid()) ambient::controller.calloc(*(revision*)this);
        return (T*)e;
    }
    template<typename T>
    inline w_revision::operator T* (){
        revision::entry& e = ((revision*)this)->content;
        if(!e.valid()){
            revision::entry& parent = this->get_parent()->content;
            if(parent.occupied()) ambient::controller.alloc(*(revision*)this);
            else e.swap(parent);
        }
        return (T*)e;
    }
    template<typename T>
    inline r_revision::operator T* (){
        revision::entry& e = ((revision*)this)->content;
        if(!e.valid()){
            revision::entry& parent = this->get_parent()->content;
            if(parent.occupied()){
                if(this->get_parent()->clean) ambient::controller.calloc(*(revision*)this);
                else ambient::controller.alloc(*(revision*)this);
            }else e.swap(parent);
        }
        return (T*)e;
    }
    
    template<class T>
    inline iteratable<T>::iteratable(dim2 dim) : T(dim) { }

    template<class T>
    inline size_t iteratable<T>::get_thread_revision_base() const {
        return this->thread_revision_base[GET_TID];
    }

    template<class T>
    inline void iteratable<T>::set_thread_revision_base(size_t r){
        this->thread_revision_base[GET_TID] = r;
    }

} } }
