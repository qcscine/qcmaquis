#ifndef AMBIENT_NUMERIC_FUTURE
#define AMBIENT_NUMERIC_FUTURE

#include <alps/numeric/real.hpp> 

#ifdef HAVE_ALPS_HDF5
#include <alps/hdf5.hpp>
#endif

namespace ambient { namespace numeric {

    using ambient::models::velvet::transformable;
    using ambient::models::velvet::transformable_expr;
    using ambient::models::velvet::transformable_value;

    template <typename T>
    class future {
    private:
        template<typename S> future& operator = (const S& v){ }
    public:
        typedef T value_type;

        void init(value_type v = T()){
            core = new (ambient::pool.malloc<FUTURE_SIZE>()) transformable_value<T>(v);
            valid = true;
        }
        template<typename S>
        void reuse(future<S>& f){
            core = (transformable*)f.core; // unsafe - proper convertion should be done
            valid = f.valid;
            f.clear();
        }
       ~future(){ if(core) ambient::destroy(core); }
        explicit constexpr future(transformable* c): core(c) {} // kernel's inner usage (no desctruction)
        template <typename FP, FP OP> explicit future(transformable_expr<T,FP,OP>* c): core(c), valid(false){
        }
        future()                                       { init();                                                            }
        future(double v)                               { init(v);                                                           }
        future(std::complex<double> v)                 { init(v);                                                           }
        future(const future& f)                        { init(f.get()); /* important */                                     }
        future(future&& f)                             { reuse(f);                                                          }
        future& operator = (const future& f)           { core->v = f.get(); return *this;                                   }
        future& operator = (future&& f)                { if(core) ambient::destroy(core); reuse(f); return *this;           }
        template<typename S> future(const future<S>& f){ init((T)f.get());                                                  }
        template<typename S> future(future<S>&& f)     { reuse(f);                                                          }
        operator T () const                            { return get();                                                      }
        const T& get_naked() const                     { return core->v;                                                    }
        T& get_naked()                                 { return core->v;                                                    }
        const future<T>& unfold() const                { assert(valid); return *this;                                       }
        future<T>& unfold()                            { assert(valid); valid = false; return *this;                        }
        T get() const                                  { if(!valid){ ambient::sync(); valid = true; } return core->eval();  }
        future& operator += (const future& r)          { valid &= r.valid; *core += *r.core; r.clear(); return *this; }
        future& operator /= (const future& r)          { core->v = get() / r.get(); return *this;                         }
        void clear() const                             { core = NULL; }
    public:
        mutable bool valid;
        mutable transformable* core;
        #ifdef HAVE_ALPS_HDF5
        void load(alps::hdf5::archive & ar){}
        void save(alps::hdf5::archive & ar) const {}
        #endif
    };

    template<typename T> future<T> operator + (const future<T>& l, const future<T>& r){
        transformable* a = l.core; l.clear();
        transformable* b = r.core; r.clear();
        return future<T>(new (ambient::pool.malloc<FUTURE_SIZE>()) 
                         transformable_expr<T, decltype(&ambient::models::velvet::op_plus<T>), 
                                            ambient::models::velvet::op_plus>(a, b)
                        ); 
    }
    #ifdef AMBIENT_LOOSE_FUTURE
    template<typename T> future<T> operator / (const future<T>& l, const future<T>& r){ 
        return future<T>(new (ambient::pool.malloc<FUTURE_SIZE>()) 
                         transformable_expr<T, decltype(&ambient::models::velvet::op_div<T>),
                                            ambient::models::velvet::op_div>(l.core, r.core)
                        ); 
    }
    inline future<double> sqrt(const future<double>& f){
        return future<double>(new (ambient::pool.malloc<FUTURE_SIZE>()) 
                              transformable_expr<double, decltype(&ambient::models::velvet::op_sqrt<double>),
                              ambient::models::velvet::op_sqrt>(f.core)
                             ); 
    }
    #else
    inline double sqrt(const future<double>& f){ return std::sqrt(f.get()); }
    template<typename T>       T operator / (const future<T>& l, const future<T>& r)    { return (l.get() / r.get());  }
    #endif
    template<typename T> const T operator / (double l, const future<T>& r)              { return (l / r.get());        }
    template<typename T> const T operator / (std::complex<double> l, const future<T>& r){ return (l / r.get());        }
    template<typename T> const future<double>& real(const future<T>& f)                 { return *(future<double>*)&f; }
    template<typename T> std::vector<double> real(const std::vector<future<T> >& f){
        ambient::sync();
        std::vector<double> res;
        for(size_t k = 0; k < f.size(); ++k) res.push_back(std::real(f[k].get()));
        return res;
    }

} }

#endif
