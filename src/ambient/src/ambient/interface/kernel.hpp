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
    class kernel_dispatch : public cfunctor
    {
    public:
        virtual ~kernel_dispatch(){
            kernel_inliner<typename K::F,K::c>::cleanup(this);
        }
        virtual void weight(){
            kernel_inliner<typename K::F,K::c>::weight(this);
        }
        virtual void place(){
            kernel_inliner<typename K::F,K::c>::place(this);
        }
        virtual void computation(){
            kernel_inliner<typename K::F,K::c>::invoke(this);
            this->check_complete();
        }
        virtual void logistics(){
            kernel_inliner<typename K::F,K::l>::invoke(this);
            this->check_complete();
        }
        virtual bool pretend(){
            return false;
        }
    };

    template<class K>
    class kernel_dispatch_unpinned : public cfunctor
    {
    public:
        virtual ~kernel_dispatch_unpinned(){
            kernel_inliner<typename K::F,K::c>::cleanup(this);
        }
        virtual void weight(){
            kernel_inliner<typename K::F,K::c>::weight(this);
        }
        virtual void place(){
            kernel_inliner<typename K::F,K::c>::place(this);
        }
        virtual void computation(){
            kernel_inliner<typename K::F,K::c>::invoke(this);
            this->check_complete();
        }
        virtual void logistics(){
            kernel_inliner<typename K::F,K::l>::invoke(this);
            this->lock();
            if(this->workload == 1) controller.execute_free_mod(this);
            else this->workload--;
            this->unlock();
        }
        virtual bool pretend(){
            this->lock();
            this->workload--;
            this->unlock();
            if(!this->workload) return false;
            return true;
        }
    };
}

#include "ambient/interface/pp/push.pp.hpp"    
#endif
