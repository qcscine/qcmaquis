#ifndef AMBIENT_INTERFACE_KERNELS
#define AMBIENT_INTERFACE_KERNELS
#include "ambient/utils/timings.hpp"

#ifdef AMBIENT_COMPUTATIONAL_TIMINGS
#include <typeinfo>
#endif

namespace ambient {

    using ambient::controllers::velvet::cfunctor;

    template<typename FP, FP fp> struct kernel_inliner{};
    #include "ambient/interface/pp/kernel_inliner.pp.hpp"

    template<class K>
    class kernel : public cfunctor {
    public:
        #define inliner kernel_inliner<decltype(&K::c),&K::c>
        inline void operator delete (void* ptr){ }
        inline void* operator new (size_t size){
            return ambient::bulk.malloc<sizeof(K)+sizeof(void*)*inliner::arity>();
        }

        kernel(){ 
            this->arguments = (void**)(this+1); // note: trashing the vtptr of derived object
            #ifdef AMBIENT_COMPUTATIONAL_DATAFLOW
            this->number = ambient::model.op_sid++;
            #endif
        }

        virtual bool ready(){ 
            return inliner::ready(this);
        }
        virtual void invoke(){ 
            #ifdef AMBIENT_COMPUTATIONAL_TIMINGS
            static ambient::timer time(typeid(K).name()); time.begin();
            #endif
            inliner::invoke(this);
            inliner::cleanup(this);
            #ifdef AMBIENT_COMPUTATIONAL_TIMINGS
            time.end();
            #endif
        }
        #ifdef AMBIENT_COMPUTATIONAL_DATAFLOW
        virtual const char* name(){ 
            return typeid(K).name(); 
        }
        #endif
        template <complexity O, class T0>
        static inline void spawn(T0& arg0){
            inliner::latch<O>(new kernel(), info<T0>::unfold<inliner::t0>(arg0));
        }
        template <complexity O, class T0, class T1>
        static inline void spawn(T0& arg0, T1& arg1){
            inliner::latch<O>(new kernel(), info<T0>::unfold<inliner::t0>(arg0) , info<T1>::unfold<inliner::t1>(arg1));
        }
        template <complexity O, class T0, class T1, class T2>
        static inline void spawn(T0& arg0, T1& arg1, T2& arg2){
            inliner::latch<O>(new kernel(), info<T0>::unfold<inliner::t0>(arg0) , info<T1>::unfold<inliner::t1>(arg1) , info<T2>::unfold<inliner::t2>(arg2));
        }
        template <complexity O, class T0 , class T1 , class T2 , class T3 >
        static inline void spawn(T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 ){
            inliner::latch<O>(new kernel(), info<T0>::unfold<inliner::t0>(arg0) , info<T1>::unfold<inliner::t1>(arg1) , info<T2>::unfold<inliner::t2>(arg2) , info<T3>::unfold<inliner::t3>(arg3) );
        }
        template <complexity O, class T0 , class T1 , class T2 , class T3 , class T4 >
        static inline void spawn(T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 , T4 &arg4 ){
            inliner::latch<O>(new kernel(), info<T0>::unfold<inliner::t0>(arg0) , info<T1>::unfold<inliner::t1>(arg1) , info<T2>::unfold<inliner::t2>(arg2) , info<T3>::unfold<inliner::t3>(arg3) , info<T4>::unfold<inliner::t4>(arg4) );
        }
        template <complexity O, class T0 , class T1 , class T2 , class T3 , class T4 , class T5 >
        static inline void spawn(T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 , T4 &arg4 , T5 &arg5 ){
            inliner::latch<O>(new kernel(), info<T0>::unfold<inliner::t0>(arg0) , info<T1>::unfold<inliner::t1>(arg1) , info<T2>::unfold<inliner::t2>(arg2) , info<T3>::unfold<inliner::t3>(arg3) , info<T4>::unfold<inliner::t4>(arg4) , info<T5>::unfold<inliner::t5>(arg5) );
        }
        template <complexity O, class T0 , class T1 , class T2 , class T3 , class T4 , class T5 , class T6 >
        static inline void spawn(T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 , T4 &arg4 , T5 &arg5 , T6 &arg6 ){
            inliner::latch<O>(new kernel(), info<T0>::unfold<inliner::t0>(arg0) , info<T1>::unfold<inliner::t1>(arg1) , info<T2>::unfold<inliner::t2>(arg2) , info<T3>::unfold<inliner::t3>(arg3) , info<T4>::unfold<inliner::t4>(arg4) , info<T5>::unfold<inliner::t5>(arg5) , info<T6>::unfold<inliner::t6>(arg6) );
        }
        template <complexity O, class T0 , class T1 , class T2 , class T3 , class T4 , class T5 , class T6 , class T7 >
        static inline void spawn(T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 , T4 &arg4 , T5 &arg5 , T6 &arg6 , T7 &arg7 ){
            inliner::latch<O>(new kernel(), info<T0>::unfold<inliner::t0>(arg0) , info<T1>::unfold<inliner::t1>(arg1) , info<T2>::unfold<inliner::t2>(arg2) , info<T3>::unfold<inliner::t3>(arg3) , info<T4>::unfold<inliner::t4>(arg4) , info<T5>::unfold<inliner::t5>(arg5) , info<T6>::unfold<inliner::t6>(arg6) , info<T7>::unfold<inliner::t7>(arg7) );
        }
        template <complexity O, class T0 , class T1 , class T2 , class T3 , class T4 , class T5 , class T6 , class T7 , class T8 >
        static inline void spawn(T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 , T4 &arg4 , T5 &arg5 , T6 &arg6 , T7 &arg7 , T8 &arg8 ){
            inliner::latch<O>(new kernel(), info<T0>::unfold<inliner::t0>(arg0) , info<T1>::unfold<inliner::t1>(arg1) , info<T2>::unfold<inliner::t2>(arg2) , info<T3>::unfold<inliner::t3>(arg3) , info<T4>::unfold<inliner::t4>(arg4) , info<T5>::unfold<inliner::t5>(arg5) , info<T6>::unfold<inliner::t6>(arg6) , info<T7>::unfold<inliner::t7>(arg7) , info<T8>::unfold<inliner::t8>(arg8) );
        }
        template <complexity O, class T0 , class T1 , class T2 , class T3 , class T4 , class T5 , class T6 , class T7 , class T8 , class T9 >
        static inline void spawn(T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 , T4 &arg4 , T5 &arg5 , T6 &arg6 , T7 &arg7 , T8 &arg8 , T9 &arg9 ){
            inliner::latch<O>(new kernel(), info<T0>::unfold<inliner::t0>(arg0) , info<T1>::unfold<inliner::t1>(arg1) , info<T2>::unfold<inliner::t2>(arg2) , info<T3>::unfold<inliner::t3>(arg3) , info<T4>::unfold<inliner::t4>(arg4) , info<T5>::unfold<inliner::t5>(arg5) , info<T6>::unfold<inliner::t6>(arg6) , info<T7>::unfold<inliner::t7>(arg7) , info<T8>::unfold<inliner::t8>(arg8) , info<T9>::unfold<inliner::t9>(arg9) );
        }
        template <complexity O, class T0 , class T1 , class T2 , class T3 , class T4 , class T5 , class T6 , class T7 , class T8 , class T9, class T10 >
        static inline void spawn(T0 &arg0 , T1 &arg1 , T2 &arg2 , T3 &arg3 , T4 &arg4 , T5 &arg5 , T6 &arg6 , T7 &arg7 , T8 &arg8 , T9 &arg9, T10 &arg10 ){
            inliner::latch<O>(new kernel(), info<T0>::unfold<inliner::t0>(arg0) , info<T1>::unfold<inliner::t1>(arg1) , info<T2>::unfold<inliner::t2>(arg2) , info<T3>::unfold<inliner::t3>(arg3) , info<T4>::unfold<inliner::t4>(arg4) , info<T5>::unfold<inliner::t5>(arg5) , info<T6>::unfold<inliner::t6>(arg6) , info<T7>::unfold<inliner::t7>(arg7) , info<T8>::unfold<inliner::t8>(arg8) , info<T9>::unfold<inliner::t9>(arg9), info<T10>::unfold<inliner::t10>(arg10) );
        }
        #undef inliner
    };
}

#endif
