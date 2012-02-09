#ifndef AMBIENT_CORE_OPERATION_H
#define AMBIENT_CORE_OPERATION_H
#include <stdlib.h>
#include <boost/shared_ptr.hpp>
#include "ambient/ambient.h"
#include "ambient/models/v_model.h"

namespace ambient { namespace models {

    // {{{ compile-time type information

    // {{{ singular types or simple types //
    template <typename T> struct singular_t_info {
        typedef T* ptr_type;
        static bool parallel(T& obj){
            return false;
        } 
        static ptr_type pointer(T& obj){
            return ptr_type(new T(obj));
        }
        static T& dereference(void* ptr){
            return *(ptr_type)ptr;
        }
        static void deallocate(void* ptr){
            delete (ptr_type)ptr;
        }
        static bool constness(T& obj){
            return false;
        }
        static void modify(T& obj, imodel::modifier* m){ 
            // empty for serial objects
        }
        static void weight(void* ptr, imodel::modifier* m){
            // empty for serial objects
        }
    };

    // types are singular by defult (see forwarding.h)
    template <typename T> struct info
    { typedef singular_t_info<T> typed; };
    // }}}

    // {{{ parallel_t derived types //
    template <typename T> struct parallel_t_info {
        typedef typename T::ptr ptr_type;
        static bool parallel(T& obj){
            return true;
        } 
        static ptr_type* pointer(T& obj){
            return new ptr_type(&obj);
        }
        static ptr_type* pointer(const T& obj){
            return new ptr_type(const_cast<T*>(&obj));
        }
        static T& dereference(void* ptr){
            return *(*(ptr_type*)ptr);
        }
        static void deallocate(void* ptr){
            delete (ptr_type*)ptr;
        }
        static void modify(const T& obj, imodel::modifier* m){
            current(obj).add_modifier(m);
        }
        static void modify(T& obj, imodel::modifier* m){
            ambient::model.add_revision(&obj);
            current(obj).add_modifier(m);
        }
        static void weight(void* ptr, imodel::modifier* m){
            T& obj = *(*(ptr_type*)ptr);
            if(obj.revisions.size() > m->get_weight()){
                m->set_weight(obj.revisions.size());
                m->set_vellum(obj);
            }
        }
        static bool constness(const T& obj){
            return true;
        }
        static bool constness(T& obj){
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
        size_t get_weight();
        void set_weight(size_t credit);
        void set_vellum(imodel::object& v);
        imodel::object& get_vellum();
    private:
        enum { MARKUP, LOGISTICS, COMPUTING } state;
        void(operation::*prototype)();
        void(operation::*cleanup)();
        void(operation::*creditup)();
        void(*logistics_ptr)();
        void(*computing_ptr)();
        void(*op)();
        void** arguments;
        size_t count;
    public:
        size_t credit;
        imodel::object* vellum;
        channels::group* grp;
        void* pin;
    };
} }

#endif
