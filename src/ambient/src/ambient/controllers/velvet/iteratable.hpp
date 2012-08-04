namespace ambient { namespace controllers { namespace velvet {

    template<typename T>
    inline c_revision::operator T* (){
        revision& e = *(revision*)this;
        if(!e.valid()){
            if(e.clean) ambient::controller.calloc(e);
            else ambient::controller.alloc(e);
        }
        return (T*)e;
    }
    template<typename T>
    inline p_revision::operator T* (){
        revision& e = *(revision*)this;
        if(!e.valid()) ambient::controller.calloc(e);
        return (T*)e;
    }
    template<typename T>
    inline w_revision::operator T* (){
        revision& e = *(revision*)this;
        if(!e.valid()){
            revision& parent = *e.get_parent();
            if(parent.occupied()) ambient::controller.alloc(e);
            else e.swap(parent);
        }
        return (T*)e;
    }
    template<typename T>
    inline r_revision::operator T* (){
        revision& e = *(revision*)this;
        if(!e.valid()){
            revision& parent = *e.get_parent();
            if(parent.occupied()){
                if(parent.clean) ambient::controller.calloc(e);
                else ambient::controller.alloc(e);
            }else e.swap(parent);
        }
        return (T*)e;
    }
    
} } }
