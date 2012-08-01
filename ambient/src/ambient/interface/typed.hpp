#ifndef AMBIENT_INTERFACE_TYPED
#define AMBIENT_INTERFACE_TYPED

namespace ambient { namespace numeric {
    template <class T> class matrix;
    template <class T> class weak_view;
    template <class T> class transpose_view;
    template <class T> class diagonal_matrix;
} }

namespace ambient { 

    using ambient::models::velvet::sfunctor;
    using ambient::controllers::velvet::cfunctor;
    using ambient::models::velvet::history;
 
    // {{{ compile-time type info: singular types or simple types
    template <typename T> struct singular_info {
        template<size_t arg> static inline void deallocate    (sfunctor* m){ } //((T*)m->arguments[arg])->~T(); 
        template<size_t arg> static inline T&   revised       (sfunctor* m){ return *(T*)m->arguments[arg]; }
        template<size_t arg> static inline void modify(T& obj, sfunctor* m){ m->arguments[arg] = (void*)new(ambient::bulk_pool.get<T>()) T(obj); }
        template<size_t arg> static inline void weight        (cfunctor* m){                                }
        template<size_t arg> static inline void place         (sfunctor* m){                                }
    };
    // }}}
    // {{{ compile-time type info: future types
    template <typename T> struct future_info {
        template<size_t arg> static inline void deallocate          (sfunctor* m){ ((T*)m->arguments[arg])->~T();            }
        template<size_t arg> static inline T&   revised             (sfunctor* m){ return *(T*)m->arguments[arg];            }
        template<size_t arg> static inline void modify(const T& obj, sfunctor* m){ m->arguments[arg] = (void*)new(ambient::bulk_pool.get<T>()) T(obj.ghost);    }
        template<size_t arg> static inline void weight              (cfunctor* m){                                           }
        template<size_t arg> static inline void place               (sfunctor* m){                                           }
    };
    // }}}
    // {{{ compile-time type info: iteratable derived types
    template <typename T> struct iteratable_info {
        template<size_t arg> 
        static inline void deallocate(sfunctor* m){
            ((T*)m->arguments[arg])->impl->clean();
            ((T*)m->arguments[arg])->~T();
        }
        template<size_t arg>
        static inline T& revised(sfunctor* m){
            return *(T*)m->arguments[arg];
        }
        template<size_t arg>
        static inline void modify(T& obj, sfunctor* m){
            history& o = *obj.impl;
            m->arguments[arg] = (void*)new(ambient::bulk_pool.get<T>()) T(obj.impl, ambient::model.time(&o));
            m->add_dependency(o.back());
            m->add_derivative(ambient::model.add_revision(&o)); 
        }
        template<size_t arg>
        static inline void modify(const T& obj, sfunctor* m){
            history& o = *obj.impl;
            m->arguments[arg] = (void*)new(ambient::bulk_pool.get<T>()) T(const_cast<T*>(&obj)->impl, ambient::model.time(&o));
            m->add_dependency(o.back());
        }
        template<size_t arg>
        static inline void weight(cfunctor* m){ }
        template<size_t arg>
        static inline void place(sfunctor* m){ }
    };
    // }}}
    // {{{ compile-time type info: weak iteratable derived types
    template <typename T> struct weak_iteratable_info {
        template<size_t arg>
        static inline void modify(T& obj, sfunctor* m){
            history& o = *obj.impl;
            m->arguments[arg] = (void*)new(ambient::bulk_pool.get<T>()) T(obj.impl, ambient::model.time(&o));
            m->add_derivative(ambient::model.add_revision(&o)); 
        }
        template<size_t arg> 
        static inline void deallocate(sfunctor* m){
            ((T*)m->arguments[arg])->impl->clean();
            ((T*)m->arguments[arg])->~T();
        }
        template<size_t arg>
        static inline T& revised(sfunctor* m){
            return *(T*)m->arguments[arg];
        }
        template<size_t arg> static inline void weight(cfunctor* m){ }
        template<size_t arg> static inline void place(sfunctor* m){ }
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
    struct info < ambient::numeric::matrix<S> > {
        typedef ambient::numeric::matrix<S> type;
        typedef iteratable_info< type > typed; 
        static inline type& unfold(type& naked){ return naked; }
        typedef S value_type;
    };

    template <typename S>
    struct info < const ambient::numeric::matrix<S> > {
        typedef ambient::numeric::matrix<S> type;
        typedef iteratable_info< type > typed; 
        static inline const type& unfold(const type& naked){ return naked; }
        typedef S value_type;
    };

    template <typename S>
    struct info < ambient::numeric::diagonal_matrix<S> > {
        typedef ambient::numeric::diagonal_matrix<S> type;
        static inline ambient::numeric::matrix<S>& unfold(type& folded){ return folded.get_data(); }
    };

    template <typename S>
    struct info < const ambient::numeric::diagonal_matrix<S> > {
        typedef const ambient::numeric::diagonal_matrix<S> type;
        static inline const ambient::numeric::matrix<S>& unfold(type& folded){ return folded.get_data(); }
    };

    template <typename S>
    struct info < ambient::numeric::weak_view<S> > {
        typedef ambient::numeric::weak_view<S> type;
        typedef weak_iteratable_info< type > typed; 
        static inline type& unfold(type& naked){ return naked; }
        typedef S value_type;
    };

    template <typename S>
    struct info < const ambient::numeric::transpose_view<ambient::numeric::matrix<S> > > {
        typedef const ambient::numeric::transpose_view<ambient::numeric::matrix<S> > type;
        static inline const ambient::numeric::matrix<S>& unfold(type& folded){ return ambient::numeric::matrix<S>(folded.impl, NULL); }
    };
    // }}}
}

#endif
