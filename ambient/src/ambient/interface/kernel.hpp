#ifndef AMBIENT_INTERFACE_KERNELS
#define AMBIENT_INTERFACE_KERNELS
#include "ambient/utils/timings.hpp"
 
//#define AMBIENT_PUSH_TIMINGS
#ifdef AMBIENT_PUSH_TIMINGS
    #define __A_TIME(name) static __a_timer time(name); time.begin();
    #define __A_TIME_STOP time.end();
#else
    #define __A_TIME(name) 
    #define __A_TIME_STOP 
#endif

namespace ambient {

    using ambient::controllers::velvet::cfunctor;
    using ambient::models::velvet::revision;

    using ambient::controllers::velvet::c_revision;
    using ambient::controllers::velvet::w_revision;
    using ambient::controllers::velvet::p_revision;
    using ambient::controllers::velvet::r_revision;

    template<typename FP, FP fp>
    struct kernel_inliner { 
        private:
        static void invoke    (){}
        static void latch     (){}
        static void cleanup   (){}
        static void place     (){}
        static bool ready(void*){}
        static bool match(void*){}
        static void tag  (void*){}
    };

    #include "ambient/interface/pp/kernel_inliner.pp.hpp"

    template<class K>
    class kernel_atomic : public cfunctor
    {
    public:
        inline void* operator new (size_t size){ 
            return ambient::bulk_pool.get<sizeof(K)>(); 
        }
        
        inline void operator delete (void* ptr){ }

        virtual void place()       { return kernel_inliner<typename K::F,&K::c>::place(this);      }
        virtual bool ready(void* e){ return kernel_inliner<typename K::F,&K::c>::ready(this, e);   }
        virtual bool match(void* t){ return kernel_inliner<typename K::F,&K::c>::match(this, t);   }
        virtual void tag(void* t)  { return kernel_inliner<typename K::F,&K::c>::tag(this, t);     }
        virtual void computation() {        kernel_inliner<typename K::F,&K::c>::invoke((K*)this); 
                                     return kernel_inliner<typename K::F,&K::c>::cleanup(this);    }
        virtual void logistics()   { return kernel_inliner<typename K::F,&K::l>::invoke((K*)this); }
        inline void ctxt_select(const char* sql){ 
            this->set_group(channel.world()); 
        }
        inline void pin(revision& r){
            ambient::controller.ifetch(r);
            ambient::controller.schedule(this);
        }
        inline void assign(revision& r){ 
            ambient::controller.ifetch(r);
        }
        template <typename T> inline revision&   ui_l_current(T& obj){ return *obj.impl->content[obj.ref];                  }
        template <typename T> inline c_revision& ui_c_current(T& obj){ return *(c_revision*)obj.impl->content[obj.ref];     }
        template <typename T> inline w_revision& ui_w_updated(T& obj){ return *(w_revision*)(obj.impl->content[obj.ref+1]); }
        template <typename T> inline p_revision& ui_p_updated(T& obj){ return *(p_revision*)(obj.impl->content[obj.ref+1]); }
        template <typename T> inline r_revision& ui_r_updated(T& obj){ return *(r_revision*)(obj.impl->content[obj.ref+1]); }
    };
}

#include "ambient/interface/pp/push.pp.hpp"    
#endif
