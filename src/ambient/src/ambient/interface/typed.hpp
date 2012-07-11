#ifndef AMBIENT_INTERFACE_TYPED
#define AMBIENT_INTERFACE_TYPED

namespace ambient { namespace numeric {
    template <class T> class matrix;
    template <class T> class matrix_impl;
    template <class T> class diagonal_matrix;
    template <class T> class transpose_view;
} }

namespace ambient { 

    using ambient::models::velvet::sfunctor;
    using ambient::controllers::velvet::cfunctor;
 
    // {{{ compile-time type info: singular types or simple types
    template <typename T> struct singular_info {
        typedef T* ptr_type;
        template<size_t arg> static inline void deallocate    (sfunctor* m){ delete  (ptr_type)(m->arguments[arg]); }
        template<size_t arg> static inline T&   revised       (sfunctor* m){ return *(ptr_type)(m->arguments[arg]); }
        template<size_t arg> static inline void modify(T& obj, sfunctor* m){ m->arguments[arg] = (void*)new T(obj); }
        template<size_t arg> static inline void weight        (cfunctor* m){                                        }
        template<size_t arg> static inline void place         (sfunctor* m){                                        }
    };
    // }}}
    // {{{ compile-time type info: future types
    template <typename T> struct future_info {
        typedef T* ptr_type;
        template<size_t arg> static inline void deallocate          (sfunctor* m){ delete (ptr_type)(m->arguments[arg]);           }
        template<size_t arg> static inline T&   revised             (sfunctor* m){ return *(ptr_type)(m->arguments[arg]);          }
        template<size_t arg> static inline void modify(const T& obj, sfunctor* m){ m->arguments[arg] = (void*)new T(obj.ghost);    }
        template<size_t arg> static inline void weight              (cfunctor* m){                                                 }
        template<size_t arg> static inline void place               (sfunctor* m){                                                 }
    };
    // }}}
    // {{{ compile-time type info: iteratable derived types
    template <typename T> struct iteratable_info {
        typedef typename T::ptr ptr_type;
        template<size_t arg> 
        static inline void deallocate(sfunctor* m){
            (*(ptr_type*)m->arguments[arg])->clean();
            delete (ptr_type*)(m->arguments[arg]);
        }
        template<size_t arg>
        static inline T& revised(sfunctor* m){
            T& obj = *(*(ptr_type*)m->arguments[arg]);
            obj.set_thread_revision_base(m->revisions[arg]);
            return obj;
        }
        template<size_t arg>
        static inline void modify(T& obj, sfunctor* m){
            m->arguments[arg] = (void*)new ptr_type(&obj);
            m->revisions[arg] = ambient::model.time(&obj);
            obj.back()->add_modifier(m);
            m->add_derivative(ambient::model.add_revision(&obj)); 
        }
        template<size_t arg>
        static inline void modify(const T& obj, sfunctor* m){
            m->arguments[arg] = (void*)new ptr_type(const_cast<T*>(&obj));
            m->revisions[arg] = ambient::model.time(&obj);
            obj.back()->add_modifier(m);
        }
        template<size_t arg>
        static inline void weight(cfunctor* m){
            T& obj = *(*(ptr_type*)(m->arguments[arg]));
            if(obj.time() > m->get_weight()){
                m->set_weight(obj.time());
                //m->set_vellum(current(obj)); // fixme
            }
        }
        template<size_t arg>
        static inline void place(sfunctor* m){
            //T& obj = *(*(ptr_type*)m->arguments[arg]);
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

    template <typename T> struct info < boost::intrusive_ptr<T> > {};
    template <typename T> struct info < const boost::intrusive_ptr<T> > {};

    template <typename S>
    struct info < ambient::future<S> > {
        typedef ambient::future<S> type;
        typedef future_info<type> typed; 
        static inline type& unfold(type& folded){ return folded.unfold(); }
    };

    template <typename S>
    struct info < const ambient::future<S> > { 
        typedef const ambient::future<S> type;
        typedef future_info<type> typed; 
        static inline type& unfold(type& folded){ return folded.unfold(); }
    };

    template <typename S>
    struct info < ambient::numeric::matrix_impl<S> > {
        typedef ambient::numeric::matrix_impl<S> type;
        typedef iteratable_info< type > typed; 
        static inline type& unfold(type& naked){ return naked; }
        typedef S value_type;
    };

    template <typename S>
    struct info < const ambient::numeric::matrix_impl<S> > {
        typedef ambient::numeric::matrix_impl<S> type;
        typedef iteratable_info< type > typed; 
        static inline const type& unfold(const type& naked){ return naked; }
        typedef S value_type;
    };

    template <typename S>
    struct info < ambient::numeric::diagonal_matrix<S> > {
        typedef ambient::numeric::diagonal_matrix<S> type;
        static inline ambient::numeric::matrix_impl<S>& unfold(type& folded){ return *folded.get_data().impl; }
    };

    template <typename S>
    struct info < const ambient::numeric::diagonal_matrix<S> > {
        typedef const ambient::numeric::diagonal_matrix<S> type;
        static inline const ambient::numeric::matrix_impl<S>& unfold(type& folded){ return *folded.get_data().impl; }
    };

    template <typename S>
    struct info < ambient::numeric::matrix<S> > {
        typedef ambient::numeric::matrix<S> type;
        static inline ambient::numeric::matrix_impl<S>& unfold(type& folded){ return *folded.impl; }
    };

    template <typename S>
    struct info < const ambient::numeric::matrix<S> > {
        typedef const ambient::numeric::matrix<S> type;
        static inline const ambient::numeric::matrix_impl<S>& unfold(type& folded){ return *folded.impl; }
    };

    template <typename S>
    struct info < const ambient::numeric::transpose_view<ambient::numeric::matrix<S> > > {
        typedef const ambient::numeric::transpose_view<ambient::numeric::matrix<S> > type;
        static inline const ambient::numeric::matrix_impl<S>& unfold(type& folded){ return *folded.impl; }
    };
    // }}}
}

#endif
