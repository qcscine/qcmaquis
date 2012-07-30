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
        static void invoke  (){}
        static void latch   (){}
        static void cleanup (){}
        static void weight  (){}
        static void place   (){}
    };

    #include "ambient/interface/pp/kernel_inliner.pp.hpp"

    template<class K>
    class kernel_atomic : public cfunctor
    {
    public:
        virtual ~kernel_atomic()  { kernel_inliner<typename K::F,&K::c>::cleanup(this);    }
        virtual void weight()     { kernel_inliner<typename K::F,&K::c>::weight(this);     }
        virtual void place()      { kernel_inliner<typename K::F,&K::c>::place(this);      }
        virtual void computation(){ kernel_inliner<typename K::F,&K::c>::invoke((K*)this); }
        virtual void logistics()  { kernel_inliner<typename K::F,&K::l>::invoke((K*)this); }
        inline void ctxt_select(const char* sql){ 
            this->set_group(channel.world()); 
        }
        inline void pin(revision& r){
            r.content.assignments.push_back(this);
            ambient::controller.ifetch(r);
        }
        inline void assign(revision& r){ 
            ambient::controller.ifetch(r);
        }
        template <size_t arg, typename T> inline revision&   ui_l_current(T& obj){ return *obj.content[this->revisions[arg]];                  }
        template <size_t arg, typename T> inline c_revision& ui_c_current(T& obj){ return *(c_revision*)obj.content[this->revisions[arg]];     }
        template <size_t arg, typename T> inline w_revision& ui_w_updated(T& obj){ return *(w_revision*)obj.content[this->revisions[arg] + 1]; }
        template <size_t arg, typename T> inline p_revision& ui_p_updated(T& obj){ return *(p_revision*)obj.content[this->revisions[arg] + 1]; }
        template <size_t arg, typename T> inline r_revision& ui_r_updated(T& obj){ return *(r_revision*)obj.content[this->revisions[arg] + 1]; }
    };
}

#include "ambient/interface/pp/push.pp.hpp"    
#endif
