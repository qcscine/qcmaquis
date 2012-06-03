#ifndef AMBIENT_INTERFACE_KERNELS
#define AMBIENT_INTERFACE_KERNELS
#include "ambient/utils/timings.hpp"
 
#define marked NULL,
#define pinned ambient::pin_marker* ,

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

    template<typename FP, FP fp>
    struct kernel_inliner { 
        private:
        static void invoke  (){}
        static void latch   (){}
        static void pin     (){}
        static void cleanup (){}
        static void weight  (){}
        static void place   (){}
    };

    class pin_marker { /* empty class for args marking */ };
    #include "ambient/interface/pp/kernel_inliner.pp.hpp"

    template< class K >
    class kernel_dispatch : public cfunctor
    {
    public:
        virtual void weight(){
            kernel_inliner<typename K::F,K::c>::weight(this);
        }
        virtual void logistics(){
            kernel_inliner<typename K::F,K::l>::invoke(this);
            kernel_inliner<typename K::F,K::l>::l_check_complete(this);
        }
        virtual void computation(){
            kernel_inliner<typename K::F,K::c>::invoke(this);
            kernel_inliner<typename K::F,K::c>::c_check_complete(this);
        }
        virtual bool pretend(){
            return kernel_inliner<typename K::F,K::c>::pretend(this);
        }
        virtual void place(){
            kernel_inliner<typename K::F,K::c>::place(this);
        }
        virtual ~kernel_dispatch(){
            kernel_inliner<typename K::F,K::c>::cleanup(this);
        }
    };
}

#include "ambient/interface/pp/push.pp.hpp"    
#undef marked
#endif
