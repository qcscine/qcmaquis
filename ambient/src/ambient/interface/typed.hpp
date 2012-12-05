#ifndef AMBIENT_INTERFACE_TYPED
#define AMBIENT_INTERFACE_TYPED

#define extract(var) T* var = (T*)m->arguments[arg];

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
    // {{{ compile-time type info: singular types + inplace and future specializations
    template <typename T> struct singular_info {
        template<size_t arg> static inline void deallocate    (sfunctor* m){                        }
        template<size_t arg> static inline void deploy        (sfunctor* m,size_t){                 }
        template<size_t arg> static inline bool pin           (cfunctor* m){ return false;          }
        template<size_t arg> static inline bool ready         (sfunctor* m){ return true;           }
        template<size_t arg> static inline T&   revised       (sfunctor* m){ extract(o); return *o; }
        template<size_t arg> static inline void modify(T& obj, sfunctor* m){ 
            m->arguments[arg] = (void*)new(ambient::bulk_pool.get<sizeof(T)>()) T(obj); 
        }
    };
    template <typename T> struct singular_inplace_info : public singular_info<T> {
        template<size_t arg> static inline T& revised(sfunctor* m){ return *(T*)&m->arguments[arg]; }
        template<size_t arg> static inline void modify(T& obj, sfunctor* m){ *(T*)&m->arguments[arg] = obj; }
    };
    template <typename T> struct future_info : public singular_info<T> {
        template<size_t arg> static inline void modify(const T& obj, sfunctor* m){ 
            m->arguments[arg] = (void*)new(ambient::bulk_pool.get<sizeof(T)>()) T(obj.ghost);
        }
    };
    // }}}
    // {{{ compile-time type info: iteratable derived types
    template <typename T> struct iteratable_info : public singular_info<T> {
        template<size_t arg> 
        static inline void deallocate(sfunctor* m){
            extract(o);
            o->impl->content[o->ref+1]->complete();
        }
        template<size_t arg>
        static inline void modify(T& obj, sfunctor* m){
            history* o = obj.impl;
            m->arguments[arg] = (void*)new(ambient::bulk_pool.get<sizeof(T)>()) T(o, ambient::model.time(o));
            ambient::model.add_revision(o, m); 
        }
        template<size_t arg> 
        static inline bool pin(cfunctor* m){ 
            extract(o);
            void* generator = o->impl->content[o->ref]->generator;
            if(generator != NULL){
                ((cfunctor*)generator)->push_back(m);
                return true;
            }
            return false;
        }
        template<size_t arg> 
        static inline bool ready(sfunctor* m){
            extract(o);
            void* generator = o->impl->content[o->ref]->generator;
            if(generator == NULL || generator == m) return true;
            return false;
        }
        template<size_t arg>
        static inline void deploy(sfunctor* m, size_t target){
            extract(o);
            ambient::controller.sync(*o->impl->content[o->ref], target);
        }
    };
    // }}}
    // {{{ compile-time type info: only read/write iteratable derived types
    template <typename T> struct read_iteratable_info : public iteratable_info<T> {
        template<size_t arg> static inline void deallocate(sfunctor* m){
            extract(o);
            o->impl->content[o->ref]->release();
        }
        template<size_t arg> static inline void modify(T& obj, sfunctor* m){
            history* o = obj.impl;
            m->arguments[arg] = (void*)new(ambient::bulk_pool.get<sizeof(T)>()) T(o, ambient::model.time(o));
            ambient::model.use_revision(o);
        }
    };
    template <typename T> struct write_iteratable_info : public iteratable_info<T> {
        template<size_t arg> static inline void deploy(sfunctor* m,size_t){        }
        template<size_t arg> static inline bool pin   (cfunctor* m){ return false; }
        template<size_t arg> static inline bool ready (sfunctor* m){ return true;  }
    };
    // }}}
    // {{{ compile-time type info: specialization for forwarded types
    template <typename T> 
    struct info { 
        typedef singular_info<T> typed; 
        static inline T& unfold(T& naked){ return naked; }
    };

    template <>
    struct info < size_t > {
        typedef size_t type;
        typedef singular_inplace_info<type> typed; 
        static inline type& unfold(type& naked){ return naked; }
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
        typedef read_iteratable_info< type > typed; 
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
        typedef write_iteratable_info< type > typed; 
        static inline type& unfold(type& naked){ return naked; }
        typedef S value_type;
    };

    template <class Matrix>
    struct info < const ambient::numeric::transpose_view<Matrix> > {
        typedef const ambient::numeric::transpose_view<Matrix> type;
        static inline const ambient::numeric::matrix<typename Matrix::value_type>& unfold(type& folded){ 
            return *(const ambient::numeric::matrix<typename Matrix::value_type>*)&folded; 
        }
    };

    template <class Matrix>
    struct info < ambient::numeric::transpose_view<Matrix> > {
        typedef ambient::numeric::transpose_view<Matrix> type;
        static inline ambient::numeric::matrix<typename Matrix::value_type>& unfold(type& folded){ 
            return *(ambient::numeric::matrix<typename Matrix::value_type>*)&folded; 
        }
    };

    // }}}
}

#undef extract
#endif
