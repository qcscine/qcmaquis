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
        template<size_t arg> static inline void modify(T& obj, sfunctor* m){ m->arguments[arg] = (void*)new(ambient::bulk_pool.get<sizeof(T)>()) T(obj); }
        template<size_t arg> static inline void place         (sfunctor* m){                                }
        template<size_t arg> static inline bool ready         (sfunctor* m, void* e){ return true;          }
        template<size_t arg> static inline bool match         (sfunctor* m, void* t){ return false;         }
        template<size_t arg> static inline void tag           (sfunctor* m, void* t){                       }
    };
    // }}}
    // {{{ compile-time type info: future types
    template <typename T> struct future_info {
        template<size_t arg> static inline void deallocate          (sfunctor* m){                                           } 
        template<size_t arg> static inline T&   revised             (sfunctor* m){ return *(T*)m->arguments[arg];            }
        template<size_t arg> static inline void modify(const T& obj, sfunctor* m){ m->arguments[arg] = (void*)new(ambient::bulk_pool.get<sizeof(T)>()) T(obj.ghost);    }
        template<size_t arg> static inline void place               (sfunctor* m){                                           }
        template<size_t arg> static inline bool ready               (sfunctor* m, void* e){ return true;                     } // might be checked
        template<size_t arg> static inline bool match               (sfunctor* m, void* t){ return false;                    }
        template<size_t arg> static inline void tag                 (sfunctor* m, void* t){                                  }
    };
    // }}}
    // {{{ compile-time type info: iteratable derived types
    template <typename T> struct iteratable_info {
        template<size_t arg> 
        static inline void deallocate(sfunctor* m){
            T* obj = (T*)m->arguments[arg];
            obj->impl->content[obj->ref+1]->reset_generator();
        }
        template<size_t arg>
        static inline T& revised(sfunctor* m){
            return *(T*)m->arguments[arg];
        }
        template<size_t arg>
        static inline void modify(T& obj, sfunctor* m){
            history& o = *obj.impl;
            m->arguments[arg] = (void*)new(ambient::bulk_pool.get<sizeof(T)>()) T(obj.impl, ambient::model.time(&o));
            ambient::model.add_revision(&o); 
        }
        template<size_t arg>
        static inline void place(sfunctor* m){ }
        template<size_t arg> 
        static inline bool ready(sfunctor* m, void* e){
            T* obj = (T*)m->arguments[arg];
            void* generator = obj->impl->content[obj->ref]->generator;
            if(generator == NULL || generator == e) return true;
            return false;
        }
        template<size_t arg> 
        static inline bool match(sfunctor* m, void* t){
            T* obj = (T*)m->arguments[arg];
            return (obj->impl->content[obj->ref]->generator == t);
        }
        template<size_t arg> 
        static inline void tag(sfunctor* m, void* t){
            T* obj = (T*)m->arguments[arg];
            obj->impl->content[obj->ref+1]->set_generator(t);
        }
    };
    // }}}
    // {{{ compile-time type info: const iteratable derived types
    template <typename T> struct const_iteratable_info {
        template<size_t arg> 
        static inline void deallocate(sfunctor* m){
        }
        template<size_t arg>
        static inline T& revised(sfunctor* m){
            return *(T*)m->arguments[arg];
        }
        template<size_t arg>
        static inline void modify(T& obj, sfunctor* m){
            history& o = *obj.impl;
            m->arguments[arg] = (void*)new(ambient::bulk_pool.get<sizeof(T)>()) T(obj.impl, ambient::model.time(&o));
        }
        template<size_t arg>
        static inline void place(sfunctor* m){ }
        template<size_t arg> 
        static inline bool ready(sfunctor* m, void* e){
            T* obj = (T*)m->arguments[arg];
            void* generator = obj->impl->content[obj->ref]->generator;
            if(generator == NULL || generator == e) return true;
            return false;
        }
        template<size_t arg> 
        static inline bool match(sfunctor* m, void* t){
            T* obj = (T*)m->arguments[arg];
            return (obj->impl->content[obj->ref]->generator == t);
        }
        template<size_t arg> 
        static inline void tag(sfunctor* m, void* t){
        }
    };
    // }}}
    // {{{ compile-time type info: weak iteratable derived types
    template <typename T> struct weak_iteratable_info {
        template<size_t arg>
        static inline void modify(T& obj, sfunctor* m){
            history& o = *obj.impl;
            m->arguments[arg] = (void*)new(ambient::bulk_pool.get<sizeof(T)>()) T(obj.impl, ambient::model.time(&o));
            ambient::model.add_revision(&o);
        }
        template<size_t arg> 
        static inline void deallocate(sfunctor* m){
            T* obj = (T*)m->arguments[arg];
            obj->impl->content[obj->ref+1]->reset_generator();
        }
        template<size_t arg>
        static inline T& revised(sfunctor* m){
            return *(T*)m->arguments[arg];
        }
        template<size_t arg> 
        static inline void tag(sfunctor* m, void* t){
            T* obj = (T*)m->arguments[arg];
            obj->impl->content[obj->ref+1]->set_generator(t);
        }
        template<size_t arg> static inline void place(sfunctor* m){ }
        template<size_t arg> static inline bool ready(sfunctor* m, void* e){ return true;  }
        template<size_t arg> static inline bool match(sfunctor* m, void* t){ return false; }
    };
    // }}}
    // {{{ compile-time type info: specialization for forwarded types
    template <typename T> 
    struct info { 
        typedef singular_info<T> typed; 
        static inline T& unfold(T& naked){ return naked; }
    };

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
        typedef const ambient::numeric::matrix<S> type;
        typedef const_iteratable_info< type > typed; 
        static inline type& unfold(type& naked){ return naked; }
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
        static inline const ambient::numeric::matrix<S>& unfold(type& folded){ return *(const ambient::numeric::matrix<S>*)&folded; }
    };
    // }}}
}

#endif
