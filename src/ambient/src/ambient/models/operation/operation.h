#ifndef AMBIENT_MODELS_OPERATION
#define AMBIENT_MODELS_OPERATION
#include <stdlib.h>
#include <boost/shared_ptr.hpp>
#include "ambient/ambient.h"
#include "ambient/models/v_model.h"
#include "ambient/utils/timings.h"

namespace ambient { namespace models {

    // {{{ compile-time type information

    template <typename T> struct info;

    // {{{ singular types or simple types //
    template <typename T> struct singular_t_info {
        typedef T* ptr_type;
        static inline bool parallel(T& obj){
            return false;
        } 
        static inline ptr_type pointer(T& obj){
            T* r = new T(obj);
            return ptr_type(r);
        }
        static inline T& dereference(void* ptr){
            return *(ptr_type)ptr;
        }
        static inline void deallocate(void* ptr){
            delete (ptr_type)ptr;
        }
        static inline bool constness(T& obj){
            return false;
        }
        static inline size_t modify(T& obj, imodel::modifier* m){ 
            return 0; // empty for serial objects
        }
        static inline void revise(void* ptr, size_t revision){
            // empty for serial objects
        }
        static inline void weight(void* ptr, imodel::modifier* m){
            // empty for serial objects
        }
        static inline void place(void* ptr, imodel::modifier* m){
            // empty for serial objects
        }
    };
    // }}}

    // {{{ parallel_t derived types //
    template <typename T> struct parallel_t_info {
        typedef typename T::ptr ptr_type;
        static inline bool parallel(T& obj){
            return true;
        } 
        static inline ptr_type* pointer(T& obj){
            return new ptr_type(&obj);
        }
        static inline ptr_type* pointer(const T& obj){
            return new ptr_type(const_cast<T*>(&obj));
        }
        static inline T& dereference(void* ptr){
            return *(*(ptr_type*)ptr);
        }
        static inline void deallocate(void* ptr){
            delete (ptr_type*)ptr;
        }
        static inline size_t modify(T& obj, imodel::modifier* m){
            size_t base = obj.get_revision_base();
            current(obj).add_modifier(m);
            ambient::model.add_revision(&obj);
            current(obj).set_generator(m); // (updated obj)
            return base;
        }
        static inline size_t modify(const T& obj, imodel::modifier* m){
            size_t base = obj.get_revision_base();
            current(obj).add_modifier(m);
            return base;
        }
        static inline void revise(void* ptr, size_t revision){
            T& obj = *(*(ptr_type*)ptr);
            ctxt.set_revision_base(&obj, revision);
        }
        static inline void weight(void* ptr, imodel::modifier* m){
            T& obj = *(*(ptr_type*)ptr);
            if(obj.revisions.size() > m->get_weight()){
                m->set_weight(obj.revisions.size());
                m->set_vellum(current(obj));
            }
        }
        static inline void place(void* ptr, imodel::modifier* m){
            T& obj = *(*(ptr_type*)ptr);
            if(current(obj).get_placement() == NULL)
                current(obj).set_placement(m->get_group());
        }
        static inline bool constness(const T& obj){
            return true;
        }
        static inline bool constness(T& obj){
            return false;
        }
    };
    // }}}

    // }}}

    class ambient_pin {}; // empty class for args marking
    class operation : public imodel::modifier {
    public:
        ~operation();
        #include "ambient/models/operation/pp/operation.pp.hpp"
        void invoke();   // executes operation (clean way)
        void weight();   // credits up the operation
        void set_group(channels::group* grp);
        channels::group* get_group();
        size_t get_weight();
        void set_weight(size_t credit);
        void set_vellum(imodel::revision& v);
        imodel::revision& get_vellum();
        imodel::revision* get_pin();
        void add_condition();
        bool pretend();
    private:
        enum { MARKUP, COMPUTING, COMPLETE } state;
        void(operation::*prototype)();
        void(operation::*cleanup)();
        void(operation::*creditup)();
        void(operation::*place)();
        void(*logistics_ptr)();
        void(*computing_ptr)();
        void(*op)();
        void** arguments;
        size_t* revisions;
        size_t count;
        long int workload; // signed for thread-safety
        pthread_mutex_t mutex;
    public:
        size_t credit;
        imodel::revision* vellum;
        imodel::revision* pin;
        channels::group* grp;
    };

} }

#endif
