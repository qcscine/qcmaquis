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
            this->affinity = r.affinity;
            r.content.assignments.push_back(this);
            ambient::controller.ifetch_block(r,0,0);
        }
        inline void assign(revision& r){ 
            ambient::controller.ifetch_block(r,0,0);
        }
    };
}

#include "ambient/interface/pp/push.pp.hpp"    
#endif
