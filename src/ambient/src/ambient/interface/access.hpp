/*
 * Ambient, License - Version 1.0 - May 3rd, 2012
 *
 * Permission is hereby granted, free of charge, to any person or organization
 * obtaining a copy of the software and accompanying documentation covered by
 * this license (the "Software") to use, reproduce, display, distribute,
 * execute, and transmit the Software, and to prepare derivative works of the
 * Software, and to permit third-parties to whom the Software is furnished to
 * do so, all subject to the following:
 *
 * The copyright notices in the Software and this entire statement, including
 * the above license grant, this restriction and the following disclaimer,
 * must be included in all copies of the Software, in whole or in part, and
 * all derivative works of the Software, unless such copies or derivative
 * works are solely in the form of machine-executable object code generated by
 * a source language processor.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

#ifndef AMBIENT_INTERFACE_ACCESS
#define AMBIENT_INTERFACE_ACCESS

#ifdef AMBIENT_PRESERVE_BULK
#define AMBIENT_RESTRICT_SPEC if(!obj.core->temporary && c.spec.region == region_t::rbulked){ c.spec.region = region_t::rstandard; }
#else
#define AMBIENT_RESTRICT_SPEC
#endif

namespace ambient {

    using ambient::models::velvet::revision;

    template <typename T> static revision& raw(T& obj){ 
        return *obj.core->content[obj.ref];     
    }

    template <typename T> static revision& safe_raw(T& obj){
        ambient::model.touch(obj.core);
        return *obj.core->current;
    }

    template <typename T> static revision& serial(T& obj){ 
        ambient::model.touch(obj.core);
        ambient::sync(); 
        revision& c = *obj.core->current;
        assert(c.state == ambient::local || c.state == ambient::common);
        if(!c.valid()){
            AMBIENT_RESTRICT_SPEC
            c.embed(T::allocator_type::calloc(c.spec));
        }
        return c;
    }

    template <typename T> static revision& current(T& obj){ 
        revision& c = *obj.core->content[obj.ref];
        if(!c.valid()){
            AMBIENT_RESTRICT_SPEC
            c.embed(T::allocator_type::calloc(c.spec));
        }
        return c;
    }

    template <typename T> static revision& updated(T& obj){ 
        revision& c = *obj.core->content[obj.ref+1]; assert(!c.valid());
        revision& p = *obj.core->content[obj.ref];
        if(p.valid() && !p.locked() && c.spec.conserves(p.spec)) c.reuse(p);
        else{
            AMBIENT_RESTRICT_SPEC
            c.embed(T::allocator_type::alloc(c.spec));
        }
        return c;
    }

    template <typename T> static revision& revised(T& obj){ 
        revision& c = *obj.core->content[obj.ref+1]; assert(!c.valid());
        revision& p = *obj.core->content[obj.ref];
        if(!p.valid()){
            AMBIENT_RESTRICT_SPEC
            c.embed(T::allocator_type::calloc(c.spec));
        }else if(!p.locked() && c.spec.conserves(p.spec)) c.reuse(p);
        else{
            AMBIENT_RESTRICT_SPEC
            c.embed(T::allocator_type::alloc(c.spec));
            memcpy((T*)c, (T*)p, p.spec.extent);
        }
        return c;
    }

    template <typename T> static revision& emptied(T& obj){
        revision& r = updated(obj);
        memset((T*)r, 0, r.spec.extent); 
        return r;
    }
}

#undef AMBIENT_RESTRICT_SPEC
#endif
