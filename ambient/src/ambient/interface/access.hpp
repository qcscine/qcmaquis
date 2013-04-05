#ifndef AMBIENT_INTERFACE_ACCESS
#define AMBIENT_INTERFACE_ACCESS

namespace ambient {

    using ambient::models::velvet::revision;

    template <typename T> static revision& raw(T& obj){ 
        return *obj.core->content[obj.ref];     
    }

    template <typename T> static revision& serial(T& obj){ 
        ambient::model.touch(obj.core);
        ambient::sync(); 
        revision& c = *obj.core->current;
        assert(ambient::model.common(c));
        if(!c.valid()) ambient::controller.calloc(c);
        return c;
    }

    template <typename T> static revision& current(T& obj){ 
        revision& c = *obj.core->content[obj.ref];
        if(!c.valid()) ambient::controller.calloc(c);
        return c;
    }

    template <typename T> static revision& updated(T& obj){ 
        revision& c = *obj.core->content[obj.ref+1]; assert(!c.valid());
        revision& p = *obj.core->content[obj.ref];
        if(!p.valid() || p.locked() || ((!p.region) && c.region))
            ambient::controller.alloc(c);
        else c.reuse(p);
        return c;
    }

    template <typename T> static revision& revised(T& obj){ 
        revision& c = *obj.core->content[obj.ref+1]; assert(!c.valid());
        revision& p = *obj.core->content[obj.ref];
        if(!p.valid()){
            ambient::controller.calloc(c);
        }else if(p.locked() || ((!p.region) && c.region)){
            ambient::controller.alloc(c);
            memcpy((T*)c, (T*)p, p.extent);
        }else{
            c.reuse(p);
        }
        return c;
    }

    template <typename T> static revision& emptied(T& obj){
        revision& r = updated(obj);
        memset((T*)r, 0, r.extent); 
        return r;
    }
}

#endif
