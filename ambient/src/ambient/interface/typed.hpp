#ifndef AMBIENT_INTERFACE_TYPED
#define AMBIENT_INTERFACE_TYPED

#define EXTRACT(var) T* var = (T*)m->arguments[arg];

namespace ambient { namespace numeric {
    template <class T, class Allocator> class matrix;
    template <class T> class transpose_view;
    template <class T> class diagonal_matrix;
} }

namespace ambient {
    template<typename T> class default_allocator;
    using ambient::controllers::velvet::cfunctor;
    using ambient::models::velvet::history;
    // {{{ compile-time type info: singular types + inplace and future specializations
    template <typename T> struct singular_info {
        template<size_t arg> static void deallocate     (cfunctor* m){                        }
        template<size_t arg> static bool pin            (cfunctor* m){ return false;          }
        template<size_t arg> static void score          (T& obj)     {                        }
        template<size_t arg> static bool ready          (cfunctor* m){ return true;           }
        template<size_t arg> static T&   revised        (cfunctor* m){ EXTRACT(o); return *o; }
        template<size_t arg> static void modify (T& obj, cfunctor* m){
            m->arguments[arg] = (void*)new(ambient::bulk.malloc<sizeof(T)>()) T(obj); 
        }
        template<size_t arg> static void modify_remote(T& obj)       {                        }
        template<size_t arg> static void modify_local(T& obj, cfunctor* m){
            m->arguments[arg] = (void*)new(ambient::bulk.malloc<sizeof(T)>()) T(obj);
        }
    };
    template <typename T> struct singular_inplace_info : public singular_info<T> {
        template<size_t arg> static T& revised(cfunctor* m){ return *(T*)&m->arguments[arg]; }
        template<size_t arg> static void modify_remote(T& obj){ }
        template<size_t arg> static void modify_local(T& obj, cfunctor* m){ *(T*)&m->arguments[arg] = obj; }
        template<size_t arg> static void modify(T& obj, cfunctor* m){ *(T*)&m->arguments[arg] = obj; }
    };
    template <typename T> struct future_info : public singular_info<T> {
        template<size_t arg> static void deallocate(cfunctor* m){       
            EXTRACT(o); o->core->generator = NULL;
        }
        template<size_t arg> static void modify_remote(T& obj){ 
            ambient::controller.rsync(obj.core);
        }
        template<size_t arg> static void modify_local(const T& obj, cfunctor* m){ 
            obj.core->generator = m;
            ambient::controller.lsync(obj.core);
            m->arguments[arg] = (void*)new(ambient::bulk.malloc<sizeof(T)>()) T(obj.core);
        }
        template<size_t arg> static void modify(const T& obj, cfunctor* m){ 
            m->arguments[arg] = (void*)new(ambient::bulk.malloc<sizeof(T)>()) T(obj.core);
        }
    };
    template <typename T> struct read_future_info : public future_info<T> {
        template<size_t arg> static void deallocate(cfunctor* m){ }
        template<size_t arg> static void modify_remote(T& obj){ }
        template<size_t arg> static void modify_local(const T& obj, cfunctor* m){
            m->arguments[arg] = (void*)new(ambient::bulk.malloc<sizeof(T)>()) T(obj.core);
        }
    };
    // }}}
    // {{{ compile-time type info: iteratable derived types
    template <typename T> struct iteratable_info : public singular_info<T> {
        template<size_t arg> 
        static void deallocate(cfunctor* m){
            EXTRACT(o);
            o->core->content[o->ref+1]->complete();
        }
        template<size_t arg>
        static void modify_remote(T& obj){
            history* o = obj.core;
            ambient::model.touch(o);
            ambient::controller.rsync(o->back());
            ambient::model.add_revision<ambient::remote>(o); 
        }
        template<size_t arg>
        static void modify_local(T& obj, cfunctor* m){
            history* o = obj.core;
            m->arguments[arg] = (void*)new(ambient::bulk.malloc<sizeof(T)>()) T(o, ambient::model.time(o));
            ambient::controller.lsync(o->back());
            ambient::model.add_revision<ambient::local>(o, m); 
        }
        template<size_t arg>
        static void modify(T& obj, cfunctor* m){
            history* o = obj.core;
            m->arguments[arg] = (void*)new(ambient::bulk.malloc<sizeof(T)>()) T(o, ambient::model.time(o));
            ambient::controller.sync(o->back());
            ambient::model.add_revision<ambient::common>(o, m); 
        }
        template<size_t arg> 
        static bool pin(cfunctor* m){ 
            EXTRACT(o);
            void* generator = o->core->content[o->ref]->generator;
            if(generator != NULL){
                ((cfunctor*)generator)->queue(m);
                return true;
            }
            return false;
        }
        template<size_t arg> 
        static void score(T& obj){
            ambient::controller.intend_fetch(obj.core);
            ambient::controller.intend_write(obj.core);
        }
        template<size_t arg> 
        static bool ready(cfunctor* m){
            EXTRACT(o);
            void* generator = o->core->content[o->ref]->generator;
            if(generator == NULL || generator == m) return true;
            return false;
        }
    };
    // }}}
    // {{{ compile-time type info: only read/write iteratable derived types
    template <typename T> struct read_iteratable_info : public iteratable_info<T> {
        template<size_t arg> static void deallocate(cfunctor* m){
            EXTRACT(o);
            o->core->content[o->ref]->release();
        }
        template<size_t arg> static void modify_remote(T& obj){
            history* o = obj.core;
            ambient::model.touch(o);
            ambient::controller.rsync(o->back());
        }
        template<size_t arg> static void modify_local(T& obj, cfunctor* m){
            history* o = obj.core;
            m->arguments[arg] = (void*)new(ambient::bulk.malloc<sizeof(T)>()) T(o, ambient::model.time(o));
            ambient::controller.lsync(o->back());
            ambient::model.use_revision(o);
        }
        template<size_t arg> static void modify(T& obj, cfunctor* m){
            history* o = obj.core;
            m->arguments[arg] = (void*)new(ambient::bulk.malloc<sizeof(T)>()) T(o, ambient::model.time(o));
            ambient::controller.sync(o->back());
            ambient::model.use_revision(o);
        }
        template<size_t arg> 
        static void score(T& obj){
            ambient::controller.intend_fetch(obj.core);
        }
    };
    template <typename T> struct write_iteratable_info : public iteratable_info<T> {
        template<size_t arg> static void modify_remote(T& obj){
            history* o = obj.core;
            ambient::model.touch(o);
            ambient::model.add_revision<ambient::remote>(o); 
        }
        template<size_t arg> static void modify_local(T& obj, cfunctor* m){
            history* o = obj.core;
            m->arguments[arg] = (void*)new(ambient::bulk.malloc<sizeof(T)>()) T(o, ambient::model.time(o));
            ambient::model.add_revision<ambient::local>(o, m); 
        }
        template<size_t arg> static void modify(T& obj, cfunctor* m){
            history* o = obj.core;
            m->arguments[arg] = (void*)new(ambient::bulk.malloc<sizeof(T)>()) T(o, ambient::model.time(o));
            ambient::model.add_revision<ambient::common>(o, m); 
        }
        template<size_t arg> static bool pin(cfunctor* m){ return false; }
        template<size_t arg> static void score(T& obj) {               
            ambient::controller.intend_write(obj.core);
        }
        template<size_t arg> static bool ready (cfunctor* m){ return true;  }
    };
    // }}}

    // {{{ compile-time type info: specialization for forwarded types
    using ambient::numeric::future;
    using ambient::numeric::matrix;
    using ambient::numeric::diagonal_matrix;
    using ambient::numeric::transpose_view;

    template <class T> struct unbound : public T {
        unbound(const typename T::ptr& p, size_t r) : T(p,r) {}
    };

    template <typename T> 
    struct info { 
        typedef singular_info<T> typed;
        template <typename U> static U& unfold(T& naked){ return naked; }
    };

    template <>
    struct info < size_t > {
        typedef size_t type;
        typedef singular_inplace_info<type> typed; 
        template <typename U> static type& unfold(type& naked){ return naked; }
    };

    template <typename S>
    struct info < future<S> > {
        typedef future<S> type;
        typedef future_info<type> typed; 
        template <typename U> static type& unfold(type& folded){ return folded.unfold(); }
    };

    template <typename S>
    struct info < const future<S> > { 
        typedef const future<S> type;
        typedef read_future_info<type> typed; 
        template <typename U> static type& unfold(type& folded){ return folded.unfold(); }
    };

    template <typename S, class A>
    struct info < matrix<S,A> > {
        typedef matrix<S,A> type;
        typedef iteratable_info< type > typed; 
        template <typename U> static U& unfold(type& naked){ return *static_cast<U*>(&naked); }
        typedef S value_type;
    };

    template <typename S, class A>
    struct info < const matrix<S,A> > {
        typedef const matrix<S,A> type;
        typedef read_iteratable_info< type > typed; 
        template <typename U> static type& unfold(type& naked){ return naked; }
        typedef S value_type;
    };

    template <typename S>
    struct info < diagonal_matrix<S> > {
        typedef diagonal_matrix<S> type;
        template <typename U> static U& unfold(type& folded){ return *static_cast<U*>(&folded.get_data()); }
    };

    template <typename S>
    struct info < const diagonal_matrix<S> > {
        typedef const diagonal_matrix<S> type;
        template <typename U> static const matrix<S, ambient::default_allocator<S> >& unfold(type& folded){ return folded.get_data(); }
    };

    template <class Matrix>
    struct info < const transpose_view<Matrix> > {
        typedef const transpose_view<Matrix> type;
        template <typename U> static const Matrix& unfold(type& folded){ return *(const Matrix*)&folded; }
    };

    template <class Matrix>
    struct info < transpose_view<Matrix> > {
        typedef transpose_view<Matrix> type;
        template <typename U> static Matrix& unfold(type& folded){ return *(Matrix*)&folded; }
    };

    template <typename T>
    struct info < unbound<T> > {
        typedef unbound<T> type;
        typedef write_iteratable_info< type > typed; 
        template <typename U> static type& unfold(type& naked){ return naked; }
        typedef typename T::value_type value_type;
    };

    // }}}
}

#undef EXTRACT
#endif
