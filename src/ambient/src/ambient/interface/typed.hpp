#ifndef AMBIENT_INTERFACE_TYPED
#define AMBIENT_INTERFACE_TYPED

namespace maquis { namespace types {
    template <class T> class p_dense_matrix;
    template <class T> class p_diagonal_matrix;
    template <class T> class p_dense_matrix_impl;
} }

namespace ambient { 

    using ambient::models::velvet::sfunctor;
    using ambient::controllers::velvet::cfunctor;
 
    // {{{ compile-time type info: singular types or simple types
    template <typename T> struct singular_info {
        typedef T* ptr_type;
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
        static inline size_t modify(T& obj, sfunctor* m){ 
            return 0; // empty for serial objects
        }
        static inline T& revised(void* ptr, size_t revisino){
            return dereference(ptr);
        }
        static inline void weight(void* ptr, cfunctor* m){
            // empty for serial objects
        }
        static inline void place(void* ptr, sfunctor* m){
            // empty for serial objects
        }
    };
    // }}}
    // {{{ compile-time type info: parallel derived types
    template <typename T> struct parallel_info {
        typedef typename T::ptr ptr_type;
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
        static inline size_t modify(T& obj, sfunctor* m){
            size_t timestamp = ambient::model.time(&obj);
            ui_m_current(obj).add_modifier(m);
            ambient::model.add_revision(&obj).set_generator(m);
            return timestamp;
        }
        static inline size_t modify(const T& obj, sfunctor* m){
            size_t timestamp = ambient::model.time(&obj);
            ui_m_current(obj).add_modifier(m);
            return timestamp;
        }
        static inline T& revised(void* ptr, size_t revision){
            T& obj = *(*(ptr_type*)ptr);
            ctxt.set_revision_base(&obj, revision);
            return dereference(ptr);
        }
        static inline void weight(void* ptr, cfunctor* m){
            T& obj = *(*(ptr_type*)ptr);
            if(obj.time() > m->get_weight()){
                m->set_weight(obj.time());
                //m->set_vellum(current(obj)); // fixme
            }
        }
        static inline void place(void* ptr, sfunctor* m){
            //T& obj = *(*(ptr_type*)ptr);
            //if(current(obj).get_placement() == NULL) // fixme
            //    current(obj).set_placement(m->get_group());
        }
    };
    // }}}
    // {{{ compile-time type info: specialization for forwarded types
    template <typename T> 
    struct info { 
        typedef singular_info<T> typed; 
        static inline T& unfold(T& naked){ return naked; }
    };

    template <typename T>
    struct info < boost::intrusive_ptr<T> > { 
    };

    template <typename T>
    struct info < const boost::intrusive_ptr<T> > { 
    };

    template <typename S>
    struct info < ambient::future<S> > { 
        typedef ambient::future<S> T;
        static inline S*& unfold(T& folded){
            return folded.unfold();
        }
    };

    template <typename S>
    struct info < const ambient::future<S> > { 
        typedef const ambient::future<S> T;
        static inline const S*& unfold(T& folded){
            assert(false); // remove if believe in const future
            return folded.unfold();
        }
    };

    template <typename S>
    struct info < maquis::types::p_diagonal_matrix<S> > { 
        typedef maquis::types::p_diagonal_matrix<S> T;
        static inline maquis::types::p_dense_matrix_impl<S>& unfold(T& folded){
            return *folded.get_data().impl;
        }
    };

    template <typename S>
    struct info < const maquis::types::p_diagonal_matrix<S> > { 
        typedef const maquis::types::p_diagonal_matrix<S> T;
        static inline const maquis::types::p_dense_matrix_impl<S>& unfold(T& folded){
            return *folded.get_data().impl;
        }
    };

    template <typename S>
    struct info < maquis::types::p_dense_matrix<S> > { 
        typedef maquis::types::p_dense_matrix<S> T;
        static inline maquis::types::p_dense_matrix_impl<S>& unfold(T& folded){
            return *folded.impl;
        }
    };

    template <typename S>
    struct info < const maquis::types::p_dense_matrix<S> > { 
        typedef const maquis::types::p_dense_matrix<S> T;
        static inline const maquis::types::p_dense_matrix_impl<S>& unfold(T& folded){
            return *folded.impl;
        }
    };

    template <typename S>
    struct info < maquis::types::p_dense_matrix_impl<S> > {
        typedef maquis::types::p_dense_matrix_impl<S> T;
        static inline T& unfold(T& naked){ return naked; }
        typedef parallel_info< T > typed; 
        typedef S value_type;
    };

    template <typename S>
    struct info < const maquis::types::p_dense_matrix_impl<S> > { 
        typedef maquis::types::p_dense_matrix_impl<S> T;
        static inline const T& unfold(const T& naked){ return naked; }
        typedef parallel_info< T > typed; 
        typedef S value_type;
    };
    // }}}
}

#endif
