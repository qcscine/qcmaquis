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

    using ambient::controllers::velvet::cfunctor;
    using ambient::models::velvet::transformable;
    using ambient::models::velvet::history;
    // {{{ compile-time type info: singular types + inplace and future specializations
    template <typename T> struct singular_info {
        template<size_t arg> static inline void deallocate     (cfunctor* m){                        }
        template<size_t arg> static inline bool pin            (cfunctor* m){ return false;          }
        template<size_t arg> static inline void score          (T& obj)     {                        }
        template<size_t arg> static inline bool ready          (cfunctor* m){ return true;           }
        template<size_t arg> static inline T&   revised        (cfunctor* m){ extract(o); return *o; }
        template<size_t arg> static inline void modify (T& obj, cfunctor* m){
            m->arguments[arg] = (void*)new(ambient::bulk.malloc<sizeof(T)>()) T(obj); 
        }
        template<size_t arg> static inline void modify_remote(T& obj)       {                        }
        template<size_t arg> static inline void modify_local(T& obj, cfunctor* m){
            m->arguments[arg] = (void*)new(ambient::bulk.malloc<sizeof(T)>()) T(obj);
        }
    };
    template <typename T> struct singular_inplace_info : public singular_info<T> {
        template<size_t arg> static inline T& revised(cfunctor* m){ return *(T*)&m->arguments[arg]; }
        template<size_t arg> static inline void modify_remote(T& obj){ }
        template<size_t arg> static inline void modify_local(T& obj, cfunctor* m){ *(T*)&m->arguments[arg] = obj; }
        template<size_t arg> static inline void modify(T& obj, cfunctor* m){ *(T*)&m->arguments[arg] = obj; }
    };
    template <typename T> struct future_info : public singular_info<T> {
        template<size_t arg> static inline void deallocate(cfunctor* m){       
            extract(o); o->core->generator = NULL;
        }
        template<size_t arg> static inline void modify_remote(T& obj){ 
            ambient::controller.rsync(obj.core);
        }
        template<size_t arg> static inline void modify_local(const T& obj, cfunctor* m){ 
            obj.core->generator = m;
            ambient::controller.lsync(obj.core);
            m->arguments[arg] = (void*)new(ambient::bulk.malloc<sizeof(T)>()) T(obj.core);
        }
        template<size_t arg> static inline void modify(const T& obj, cfunctor* m){ 
            m->arguments[arg] = (void*)new(ambient::bulk.malloc<sizeof(T)>()) T(obj.core);
        }
    };
    template <typename T> struct read_future_info : public future_info<T> {
        template<size_t arg> static inline void deallocate(cfunctor* m){ }
        template<size_t arg> static inline void modify_remote(T& obj){ }
        template<size_t arg> static inline void modify_local(const T& obj, cfunctor* m){
            m->arguments[arg] = (void*)new(ambient::bulk.malloc<sizeof(T)>()) T(obj.core);
        }
    };
    // }}}
    // {{{ compile-time type info: iteratable derived types
    template <typename T> struct iteratable_info : public singular_info<T> {
        template<size_t arg> 
        static inline void deallocate(cfunctor* m){
            extract(o);
            o->core->content[o->ref+1]->complete();
        }
        template<size_t arg>
        static inline void modify_remote(T& obj){
            history* o = obj.core;
            ambient::model.touch(o);
            ambient::controller.rsync(o->back());
            ambient::model.add_revision<ambient::remote>(o); 
        }
        template<size_t arg>
        static inline void modify_local(T& obj, cfunctor* m){
            history* o = obj.core;
            m->arguments[arg] = (void*)new(ambient::bulk.malloc<sizeof(T)>()) T(o, ambient::model.time(o));
            ambient::controller.lsync(o->back());
            ambient::model.add_revision<ambient::local>(o, m); 
        }
        template<size_t arg>
        static inline void modify(T& obj, cfunctor* m){
            history* o = obj.core;
            m->arguments[arg] = (void*)new(ambient::bulk.malloc<sizeof(T)>()) T(o, ambient::model.time(o));
            ambient::controller.sync(o->back());
            ambient::model.add_revision<ambient::common>(o, m); 
        }
        template<size_t arg> 
        static inline bool pin(cfunctor* m){ 
            extract(o);
            void* generator = o->core->content[o->ref]->generator;
            if(generator != NULL){
                ((cfunctor*)generator)->queue(m);
                return true;
            }
            return false;
        }
        template<size_t arg> 
        static inline void score(T& obj){
            ambient::controller.intend_fetch(obj.core);
            ambient::controller.intend_write(obj.core);
        }
        template<size_t arg> 
        static inline bool ready(cfunctor* m){
            extract(o);
            void* generator = o->core->content[o->ref]->generator;
            if(generator == NULL || generator == m) return true;
            return false;
        }
    };
    // }}}
    // {{{ compile-time type info: only read/write iteratable derived types
    template <typename T> struct read_iteratable_info : public iteratable_info<T> {
        template<size_t arg> static inline void deallocate(cfunctor* m){
            extract(o);
            o->core->content[o->ref]->release();
        }
        template<size_t arg> static inline void modify_remote(T& obj){
            history* o = obj.core;
            ambient::model.touch(o);
            ambient::controller.rsync(o->back());
        }
        template<size_t arg> static inline void modify_local(T& obj, cfunctor* m){
            history* o = obj.core;
            m->arguments[arg] = (void*)new(ambient::bulk.malloc<sizeof(T)>()) T(o, ambient::model.time(o));
            ambient::controller.lsync(o->back());
            ambient::model.use_revision(o);
        }
        template<size_t arg> static inline void modify(T& obj, cfunctor* m){
            history* o = obj.core;
            m->arguments[arg] = (void*)new(ambient::bulk.malloc<sizeof(T)>()) T(o, ambient::model.time(o));
            ambient::controller.sync(o->back());
            ambient::model.use_revision(o);
        }
        template<size_t arg> 
        static inline void score(T& obj){
            ambient::controller.intend_fetch(obj.core);
        }
    };
    template <typename T> struct write_iteratable_info : public iteratable_info<T> {
        template<size_t arg> static inline void modify_remote(T& obj){
            history* o = obj.core;
            ambient::model.touch(o);
            ambient::model.add_revision<ambient::remote>(o); 
        }
        template<size_t arg> static inline void modify_local(T& obj, cfunctor* m){
            history* o = obj.core;
            m->arguments[arg] = (void*)new(ambient::bulk.malloc<sizeof(T)>()) T(o, ambient::model.time(o));
            ambient::model.add_revision<ambient::local>(o, m); 
        }
        template<size_t arg> static inline void modify(T& obj, cfunctor* m){
            history* o = obj.core;
            m->arguments[arg] = (void*)new(ambient::bulk.malloc<sizeof(T)>()) T(o, ambient::model.time(o));
            ambient::model.add_revision<ambient::common>(o, m); 
        }
        template<size_t arg> static inline bool pin(cfunctor* m){ return false; }
        template<size_t arg> static inline void score(T& obj) {               
            ambient::controller.intend_write(obj.core);
        }
        template<size_t arg> static inline bool ready (cfunctor* m){ return true;  }
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
    struct info < ambient::numeric::future<S> > {
        typedef ambient::numeric::future<S> type;
        typedef future_info<type> typed; 
        static inline type& unfold(type& folded){ return folded.unfold(); }
    };

    template <typename S>
    struct info < const ambient::numeric::future<S> > { 
        typedef const ambient::numeric::future<S> type;
        typedef read_future_info<type> typed; 
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
