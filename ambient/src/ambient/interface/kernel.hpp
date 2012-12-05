#ifndef AMBIENT_INTERFACE_KERNELS
#define AMBIENT_INTERFACE_KERNELS
#include "ambient/utils/timings.hpp"
 
//#define AMBIENT_SPAWN_TIMINGS
#ifdef AMBIENT_SPAWN_TIMINGS
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
    using ambient::controllers::velvet::s_revision;
    using ambient::controllers::velvet::w_revision;
    using ambient::controllers::velvet::p_revision;

    template<typename FP, FP fp> struct kernel_inliner{};
    #include "ambient/interface/pp/kernel_inliner.pp.hpp"

    static inline void assign(revision& r){ } //ambient::controller.ifetch(r); }
    template <typename T> static inline revision&     current(T& obj){ return *(revision  *)obj.impl->content[obj.ref];   }
    template <typename T> static inline w_revision&   updated(T& obj){ return *(w_revision*)obj.impl->content[obj.ref+1]; }
    // supplementary revision modes: checked current (calloc), same updated (memcpy), purged updated (memset)
    template <typename T> static inline c_revision& c_current(T& obj){ return *(c_revision*)obj.impl->content[obj.ref];   }
    template <typename T> static inline s_revision& s_updated(T& obj){ return *(s_revision*)obj.impl->content[obj.ref+1]; }
    template <typename T> static inline p_revision& p_updated(T& obj){ return *(p_revision*)obj.impl->content[obj.ref+1]; }

    template<class K>
    class kernel : public cfunctor
    {
    public:
        inline void* operator new (size_t size){ 
            return ambient::bulk_pool.get<sizeof(K)>(); 
        }
        
        inline void operator delete (void* ptr){ }

        virtual bool ready()       { return kernel_inliner<typename K::F,&K::c>::ready(this);      }
        virtual bool match(void* t){ return kernel_inliner<typename K::F,&K::c>::match(this, t);   }
        virtual void place()       {        kernel_inliner<typename K::F,&K::c>::place(this);      }
        virtual void invoke()      {        kernel_inliner<typename K::F,&K::c>::invoke((K*)this); 
                                            kernel_inliner<typename K::F,&K::c>::cleanup(this);    }

        inline void pin(revision& r){
            //this->set_group(channel.world()); 
            //ambient::controller.ifetch(r);
            cfunctor* g = (cfunctor*)r.get_generator();
            if(g != NULL) g->push_back(this);
            else ambient::controller.submit(this);
        }

        template <class T0>
        static inline void spawn(T0& arg0){
            K* op = new K();
            kernel_inliner<typename K::F,&K::c>::latch(op, info<T0>::unfold(arg0));
            kernel_inliner<typename K::F,&K::l>::invoke(op); 
        }
        template <class T0, class T1>
        static inline void spawn(T0& arg0, T1& arg1){
            K* op = new K();
            kernel_inliner<typename K::F,&K::c>::latch(op, info<T0>::unfold(arg0), info<T1>::unfold(arg1));
            kernel_inliner<typename K::F,&K::l>::invoke(op); 
        }
        template <class T0, class T1, class T2>
        static inline void spawn(T0& arg0, T1& arg1, T2& arg2){
            K* op = new K();
            kernel_inliner<typename K::F,&K::c>::latch(op, info<T0>::unfold(arg0), info<T1>::unfold(arg1), info<T2>::unfold(arg2));
            kernel_inliner<typename K::F,&K::l>::invoke(op); 
        }
        template <class T0 , class T1 , class T2 , class T3 >
        static inline void spawn(T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 ){
            K* op = new K();
            kernel_inliner<typename K::F, &K::c>::latch(op, info<T0>::unfold(arg0) , info<T1>::unfold(arg1) , info<T2>::unfold(arg2) , info<T3>::unfold(arg3) );
            kernel_inliner<typename K::F,&K::l>::invoke(op); 
        }
        template <class T0 , class T1 , class T2 , class T3 , class T4 >
        static inline void spawn(T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 , T4 &arg4 ){
            K* op = new K();
            kernel_inliner<typename K::F,&K::c>::latch(op, info<T0>::unfold(arg0) , info<T1>::unfold(arg1) , info<T2>::unfold(arg2) , info<T3>::unfold(arg3) , info<T4>::unfold(arg4) );
            kernel_inliner<typename K::F,&K::l>::invoke(op); 
        }
        template <class T0 , class T1 , class T2 , class T3 , class T4 , class T5 >
        static inline void spawn(T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 , T4 &arg4 , T5 &arg5 ){
            K* op = new K();
            kernel_inliner<typename K::F,&K::c>::latch(op, info<T0>::unfold(arg0) , info<T1>::unfold(arg1) , info<T2>::unfold(arg2) , info<T3>::unfold(arg3) , info<T4>::unfold(arg4) , info<T5>::unfold(arg5) );
            kernel_inliner<typename K::F,&K::l>::invoke(op); 
        }
        template <class T0 , class T1 , class T2 , class T3 , class T4 , class T5 , class T6 >
        static inline void spawn(T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 , T4 &arg4 , T5 &arg5 , T6 &arg6 ){
            K* op = new K();
            kernel_inliner<typename K::F,&K::c>::latch(op, info<T0>::unfold(arg0) , info<T1>::unfold(arg1) , info<T2>::unfold(arg2) , info<T3>::unfold(arg3) , info<T4>::unfold(arg4) , info<T5>::unfold(arg5) , info<T6>::unfold(arg6) );
            kernel_inliner<typename K::F,&K::l>::invoke(op); 
        }
        template <class T0 , class T1 , class T2 , class T3 , class T4 , class T5 , class T6 , class T7 >
        static inline void spawn(T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 , T4 &arg4 , T5 &arg5 , T6 &arg6 , T7 &arg7 ){
            K* op = new K();
            kernel_inliner<typename K::F,&K::c>::latch(op, info<T0>::unfold(arg0) , info<T1>::unfold(arg1) , info<T2>::unfold(arg2) , info<T3>::unfold(arg3) , info<T4>::unfold(arg4) , info<T5>::unfold(arg5) , info<T6>::unfold(arg6) , info<T7>::unfold(arg7) );
            kernel_inliner<typename K::F,&K::l>::invoke(op); 
        }
        template <class T0 , class T1 , class T2 , class T3 , class T4 , class T5 , class T6 , class T7 , class T8 >
        static inline void spawn(T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 , T4 &arg4 , T5 &arg5 , T6 &arg6 , T7 &arg7 , T8 &arg8 ){
            K* op = new K();
            kernel_inliner<typename K::F,&K::c>::latch(op, info<T0>::unfold(arg0) , info<T1>::unfold(arg1) , info<T2>::unfold(arg2) , info<T3>::unfold(arg3) , info<T4>::unfold(arg4) , info<T5>::unfold(arg5) , info<T6>::unfold(arg6) , info<T7>::unfold(arg7) , info<T8>::unfold(arg8) );
            kernel_inliner<typename K::F,&K::l>::invoke(op); 
        }
        template <class T0 , class T1 , class T2 , class T3 , class T4 , class T5 , class T6 , class T7 , class T8 , class T9 >
        static inline void spawn(T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 , T4 &arg4 , T5 &arg5 , T6 &arg6 , T7 &arg7 , T8 &arg8 , T9 &arg9 ){
            K* op = new K();
            kernel_inliner<typename K::F,&K::c>::latch(op, info<T0>::unfold(arg0) , info<T1>::unfold(arg1) , info<T2>::unfold(arg2) , info<T3>::unfold(arg3) , info<T4>::unfold(arg4) , info<T5>::unfold(arg5) , info<T6>::unfold(arg6) , info<T7>::unfold(arg7) , info<T8>::unfold(arg8) , info<T9>::unfold(arg9) );
            kernel_inliner<typename K::F,&K::l>::invoke(op); 
        }
        template <class T0 , class T1 , class T2 , class T3 , class T4 , class T5 , class T6 , class T7 , class T8 , class T9, class T10 >
        static inline void spawn(T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 , T4 &arg4 , T5 &arg5 , T6 &arg6 , T7 &arg7 , T8 &arg8 , T9 &arg9, T10 &arg10 ){
            K* op = new K();
            kernel_inliner<typename K::F,&K::c>::latch(op, info<T0>::unfold(arg0) , info<T1>::unfold(arg1) , info<T2>::unfold(arg2) , info<T3>::unfold(arg3) , info<T4>::unfold(arg4) , info<T5>::unfold(arg5) , info<T6>::unfold(arg6) , info<T7>::unfold(arg7) , info<T8>::unfold(arg8) , info<T9>::unfold(arg9), info<T10>::unfold(arg10) );
            kernel_inliner<typename K::F,&K::l>::invoke(op); 
        }


    };
}

#endif
