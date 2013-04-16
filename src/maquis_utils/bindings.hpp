#ifndef MATRIX_BINDINGS_H
#define MATRIX_BINDINGS_H

#include <alps/numeric/real.hpp> 

namespace maquis { namespace bindings {

    template <typename O, typename I> struct binding { 
        static O convert(const I& m){ return static_cast<O>(m); }
    };

    template<typename O, typename I> O matrix_cast(I const& input){
       return binding<O,I>::convert(input);
    }

    template <typename T>
    struct binding< std::vector<T>, alps::numeric::diagonal_matrix<T> > {
        static std::vector<T> convert(const alps::numeric::diagonal_matrix<T>& m){
            return m.get_values();
        }
    };

    template <typename T, typename S, template<class M, class SS> class C>
    struct binding< std::vector< std::vector<T> >, C<alps::numeric::diagonal_matrix<T>, S> > {
        static std::vector< std::vector<T> > convert(const C<alps::numeric::diagonal_matrix<T>, S>& m){
            std::vector< std::vector<T> > set;
            for(size_t k = 0; k < m.n_blocks(); ++k){
                set.push_back(m[k].get_values());
            }
            return set;
        }
    };

#ifdef AMBIENT 

    using ambient::complexity;

    template <typename T, typename D>
    struct binding< ambient::numeric::tiles<ambient::numeric::diagonal_matrix<T> >, ambient::numeric::tiles<ambient::numeric::diagonal_matrix<D> > >{
       static  ambient::numeric::tiles<ambient::numeric::diagonal_matrix<T> >convert(const ambient::numeric::tiles<ambient::numeric::diagonal_matrix<D> >& m){
           ambient::numeric::tiles<ambient::numeric::diagonal_matrix<T> > pm(num_rows(m), num_cols(m));
           for(size_t k = 0; k < m.data.size(); ++k)
               ambient::numeric::kernels::cast_double_complex<T,D>::spawn<complexity::N2>(pm[k],m[k]);
           return pm;
       } 
    };

    template <typename T, typename S, template<class M, class SS> class C>
    struct binding< std::vector< std::vector<T> >, C<ambient::numeric::diagonal_matrix<T>, S> > {
        static std::vector< std::vector<T> > convert(const C<ambient::numeric::diagonal_matrix<T>, S>& m){
            std::vector< std::vector<T> > set;
            for(size_t k = 0; k < m.n_blocks(); ++k) 
                set.push_back(std::vector<T>(m[k].num_rows()));
            size_t num_cols(1);
            size_t offset(0);
            size_t num_rows;
            for(size_t k = 0; k < m.n_blocks(); ++k){
                num_rows = m[k].num_rows();
                std::vector<T>* v_ptr = &set[k];
                ambient::numeric::kernels::cast_to_vector<T>::spawn<complexity::N2>(v_ptr, m[k], num_rows, num_cols, num_rows, offset);
            }
            ambient::sync();
            return set;
        }
    };

    template <typename T>
    struct binding< ambient::numeric::matrix<T>, alps::numeric::matrix<T> > {
        static ambient::numeric::matrix<T> convert(const alps::numeric::matrix<T>& m){
            size_t num_rows = m.num_rows();
            size_t num_cols = m.num_cols();
            size_t lda = m.stride2();
            size_t offset(0);
            ambient::numeric::matrix<T> pm(num_rows, num_cols);    
            const std::vector<typename alps::numeric::matrix<T>::value_type>* v_ptr = &m.get_values();
            ambient::numeric::kernels::cast_from_vector<T>::spawn<complexity::N2>(v_ptr, pm, num_rows, num_cols, lda, offset);
            ambient::sync();
            return pm;
        }
    };

    template <typename T>
    struct binding< alps::numeric::matrix<T>, ambient::numeric::matrix<T> > {
        static alps::numeric::matrix<T> convert(const ambient::numeric::matrix<T>& pm){
            size_t num_rows = pm.num_rows();
            size_t num_cols = pm.num_cols();
            size_t offset(0);
            alps::numeric::matrix<T> m(num_rows, num_cols);    
            std::vector<typename alps::numeric::matrix<T>::value_type>* v_ptr = &m.get_values();
            ambient::numeric::kernels::cast_to_vector<T>::spawn<complexity::N2>(v_ptr, pm, num_rows, num_cols, num_rows, offset);
            ambient::sync();
            return m;
        }
    };

    template <typename T>
    struct binding< ambient::numeric::diagonal_matrix<T>, alps::numeric::diagonal_matrix<T> > {
        static ambient::numeric::diagonal_matrix<T> convert(const alps::numeric::diagonal_matrix<T>& m){
            size_t num_rows(m.num_rows());
            size_t num_cols(1);
            size_t offset(0);
            ambient::numeric::diagonal_matrix<T> pm(num_rows);    
            const std::vector<typename alps::numeric::diagonal_matrix<T>::value_type>* v_ptr = &m.get_values();
            ambient::numeric::kernels::cast_from_vector<T>::spawn<complexity::N2>(v_ptr, pm, num_rows, num_cols, num_rows, offset);
            ambient::sync();
            return pm;
        }
    };

    template <typename T>
    struct binding< alps::numeric::diagonal_matrix<T>, ambient::numeric::diagonal_matrix<T> > {
        static alps::numeric::diagonal_matrix<T> convert(const ambient::numeric::diagonal_matrix<T>& pm){
            size_t offset(0);
            size_t num_cols(1);
            size_t num_rows = pm.num_rows();
            alps::numeric::diagonal_matrix<T> m((std::size_t)num_rows);
            std::vector<typename alps::numeric::diagonal_matrix<T>::value_type>* v_ptr = &m.get_values();
            ambient::numeric::kernels::cast_to_vector<T>::spawn<complexity::N2>(v_ptr, pm, num_rows, num_cols, num_rows, offset);
            ambient::sync();
            return m;
        }
    };

    template <typename T>
    struct binding< std::vector<T>, ambient::numeric::diagonal_matrix<T> > {
        static std::vector<T> convert(const ambient::numeric::diagonal_matrix<T>& pm){
            return binding<alps::numeric::diagonal_matrix<T>, ambient::numeric::diagonal_matrix<T> >::convert(pm).get_values();
        }
    };

    template <typename T, typename S, template<class M, class SS> class C>
    struct binding< std::vector< std::vector<T> >, C<ambient::numeric::tiles<ambient::numeric::diagonal_matrix<T> >, S> > {
        static std::vector< std::vector<T> > convert(const C<ambient::numeric::tiles<ambient::numeric::diagonal_matrix<T> >, S>& m){
            std::vector< std::vector<T> > set;
            for(size_t k = 0; k < m.n_blocks(); ++k) 
                set.push_back(std::vector<T>(m[k].num_rows()));
            size_t num_cols(1);
            for(size_t k = 0; k < m.n_blocks(); ++k){
                std::vector<T>* v_ptr = &set[k];
                size_t offset = 0;
                for(size_t kk = 0; kk < m[k].data.size(); kk++){
                    size_t num_rows = m[k][kk].num_rows();
                    ambient::numeric::kernels::cast_to_vector<T>::spawn<complexity::N2>(v_ptr, m[k][kk], num_rows, num_cols, num_rows, offset);
                    offset += num_rows;
                }
            }
            ambient::sync();
            return set;
        }
    };

    template <typename T>
    struct binding< ambient::numeric::tiles<ambient::numeric::matrix<T> >, alps::numeric::matrix<T> > {
        static ambient::numeric::tiles<ambient::numeric::matrix<T> > convert(const alps::numeric::matrix<T>& m){
            size_t num_rows = m.num_rows();
            size_t num_cols = m.num_cols();
            size_t lda = m.stride2();
            ambient::numeric::tiles<ambient::numeric::matrix<T> > pm(num_rows, num_cols);    
            const std::vector<typename alps::numeric::matrix<T>::value_type>* v_ptr = &m.get_values();

            for(size_t j = 0; j < pm.nt; ++j){
                size_t offset = j*lda*AMBIENT_IB;
                for(size_t i = 0; i < pm.mt; ++i){
                    ambient::numeric::matrix<T>& tile = pm.tile(i,j);
                    size_t rows = tile.num_rows();
                    size_t cols = tile.num_cols();
                    ambient::numeric::kernels::cast_from_vector<T>::spawn<complexity::N2>(v_ptr, tile, rows, cols, lda, offset);
                    offset += rows;
                }
            }

            ambient::sync();
            return pm;
        }
    };

    template <typename T2, typename T1>
    struct binding< ambient::numeric::tiles<ambient::numeric::matrix<T2> >, alps::numeric::matrix<T1> > {
        static ambient::numeric::tiles<ambient::numeric::matrix<T2> > convert(const alps::numeric::matrix<T1>& m){
            size_t num_rows = m.num_rows();
            size_t num_cols = m.num_cols();
            size_t lda = m.stride2();
            ambient::numeric::tiles<ambient::numeric::matrix<T2> > pm(num_rows, num_cols);    
            const std::vector<typename alps::numeric::matrix<T1>::value_type>* v_ptr = &m.get_values();

            for(size_t j = 0; j < pm.nt; ++j){
                size_t offset = j*lda*AMBIENT_IB;
                for(size_t i = 0; i < pm.mt; ++i){
                    ambient::numeric::matrix<T2>& tile = pm.tile(i,j);
                    size_t rows = tile.num_rows();
                    size_t cols = tile.num_cols();
                    ambient::numeric::kernels::cast_from_vector_t<T1,T2>::spawn<complexity::N2>(v_ptr, tile, rows, cols, lda, offset);
                    offset += rows;
                }
            }

            ambient::sync();
            return pm;
        }
    };

    template <typename T>
    struct binding< alps::numeric::matrix<T>, ambient::numeric::tiles<ambient::numeric::matrix<T> > > {
        static alps::numeric::matrix<T> convert(const ambient::numeric::tiles<ambient::numeric::matrix<T> >& pm){
            size_t num_rows = pm.num_rows();
            size_t num_cols = pm.num_cols();
            alps::numeric::matrix<T> m(num_rows, num_cols);    
            std::vector<typename alps::numeric::matrix<T>::value_type>* v_ptr = &m.get_values();
            size_t lda = m.stride2();
            for(size_t j = 0; j < pm.nt; ++j){
                size_t offset = j*lda*AMBIENT_IB;
                for(size_t i = 0; i < pm.mt; ++i){
                    const ambient::numeric::matrix<T>& tile = pm.tile(i,j);
                    size_t rows = tile.num_rows();
                    size_t cols = tile.num_cols();
                    ambient::numeric::kernels::cast_to_vector<T>::spawn<complexity::N2>(v_ptr, tile, rows, cols, lda, offset);
                    offset += rows;
                }
            }

            ambient::sync();
            return m;
        }
    };

    template <typename T>
    struct binding< ambient::numeric::tiles<ambient::numeric::diagonal_matrix<T> >, alps::numeric::diagonal_matrix<T> > {
        static ambient::numeric::tiles<ambient::numeric::diagonal_matrix<T> > convert(const alps::numeric::diagonal_matrix<T>& m){
            size_t num_rows = m.num_rows();
            size_t num_cols(1);
            ambient::numeric::tiles<ambient::numeric::diagonal_matrix<T> > pm(num_rows);    
            const std::vector<typename alps::numeric::diagonal_matrix<T>::value_type>* v_ptr = &m.get_values();

            size_t offset(0);
            for(size_t i = 0; i < pm.nt; ++i){
                ambient::numeric::diagonal_matrix<T>& tile = pm[i];
                size_t rows = tile.num_rows();
                ambient::numeric::kernels::cast_from_vector<T>::spawn<complexity::N2>(v_ptr, tile, rows, num_cols, num_rows, offset);
                offset += rows;
            }

            ambient::sync();
            return pm;
        }
    };

    template <typename T>
    struct binding< alps::numeric::diagonal_matrix<T>, ambient::numeric::tiles<ambient::numeric::diagonal_matrix<T> > > {
        static alps::numeric::diagonal_matrix<T> convert(const ambient::numeric::tiles<ambient::numeric::diagonal_matrix<T> >& pm){
            size_t num_rows = pm.num_rows();
            size_t num_cols(1);
            size_t offset(0);
            alps::numeric::diagonal_matrix<T> m((std::size_t)num_rows);    
            std::vector<typename alps::numeric::diagonal_matrix<T>::value_type>* v_ptr = &m.get_values();
            for(size_t i = 0; i < pm.nt; ++i){
                const ambient::numeric::diagonal_matrix<T>& tile = pm[i];
                size_t rows = tile.num_rows();
                ambient::numeric::kernels::cast_to_vector<T>::spawn<complexity::N2>(v_ptr, tile, rows, num_cols, num_rows, offset);
                offset += rows;
            }

            ambient::sync();
            return m;
        }
    };

#endif

} }

#ifdef AMBIENT

template<typename T>
bool operator == (alps::numeric::matrix<T> const & m, ambient::numeric::matrix<T> const & pm){
    ambient::numeric::matrix<T> pm_ = maquis::bindings::matrix_cast<ambient::numeric::matrix<T> >(m);
    return (pm_ == pm);
}

template<typename T>
bool operator == (alps::numeric::diagonal_matrix<T> const & m, ambient::numeric::diagonal_matrix<T> const & pm){
    ambient::numeric::diagonal_matrix<T> pm_ = maquis::bindings::matrix_cast<ambient::numeric::diagonal_matrix<T> >(m);
    return (pm_ == pm);
}

template<typename T>
bool operator == (alps::numeric::matrix<T> const & m, ambient::numeric::tiles<ambient::numeric::matrix<T> > const & pm){
    ambient::numeric::tiles<ambient::numeric::matrix<T> > pm_ = maquis::bindings::matrix_cast<ambient::numeric::tiles<ambient::numeric::matrix<T> > >(m);
    return (pm_ == pm);
}

template<typename T>
bool operator == (alps::numeric::diagonal_matrix<T> const & m, ambient::numeric::tiles<ambient::numeric::diagonal_matrix<T> > const & pm){
    ambient::numeric::tiles<ambient::numeric::diagonal_matrix<T> > pm_ = maquis::bindings::matrix_cast<ambient::numeric::tiles<ambient::numeric::diagonal_matrix<T> > >(m);
    return (pm_ == pm);
}

template<typename T>
bool operator == (ambient::numeric::matrix<T> const & pm, alps::numeric::matrix<T> const & m){
    return (m == pm);
}

template<typename T>
bool operator == (ambient::numeric::diagonal_matrix<T> const & pm, alps::numeric::diagonal_matrix<T> const & m){
    return (m == pm);
}

template<typename T>
bool operator == (ambient::numeric::tiles<ambient::numeric::matrix<T> > const & pm, alps::numeric::matrix<T> const & m){
    return (m == pm);
}

template<typename T>
bool operator == (ambient::numeric::tiles<ambient::numeric::diagonal_matrix<T> > const & pm, alps::numeric::diagonal_matrix<T> const & m){
    return (m == pm);
}

#endif

#ifdef AMBIENT
namespace maquis {

    template<class T>
    struct real_type {
        typedef T type;
    };

    template<>
    struct real_type<std::complex<double> > {
        typedef double type;
    };

    template<>
    struct real_type<ambient::numeric::future<std::complex<double> > > {
        typedef typename ambient::numeric::future<double> type;
    };

    inline double real(double f){
        return f;
    }

    inline double real(const std::complex<double>& f){
        return f.real();
    }

    template<typename T>
    inline const typename real_type<ambient::numeric::future<T> >::type& 
    real(const ambient::numeric::future<T>& f){
        return ambient::numeric::real(f);
    }

    template <typename T,                          // value_type
             template<class AT> class A,           // allocator 
             template<class TT, class AA> class C> // vector<value_type, allocator>
    inline const C<typename real_type<T>::type, A<typename real_type<T>::type> >& 
    real(const C<ambient::numeric::future<T>, A<ambient::numeric::future<T> > >& f){
        return ambient::numeric::real(f);
    }

    template<typename _InputIterator, typename _Tp>
    inline _Tp
    accumulate(_InputIterator __first, _InputIterator __last, _Tp __init)
    {
      __glibcxx_function_requires(_InputIteratorConcept<_InputIterator>)
      __glibcxx_requires_valid_range(__first, __last);

      for (; __first != __last; ++__first)
	__init = __init + *__first;

      #ifdef AMBIENT_LOOSE_FUTURE
      return std::move(__init);
      #else
      return __init;
      #endif
    }

    template<typename T>
    inline T sqrt(T arg){
        return std::sqrt(arg);
    }

    inline ambient::numeric::future<double> sqrt(const ambient::numeric::future<double>& f){
        return ambient::numeric::sqrt(f);
    }

}
namespace ietl {

    inline double real(const ambient::numeric::future<std::complex<double> >& f){
        return ambient::numeric::real(f);
    }

}

#else

namespace maquis {

    template <class T> 
    inline typename alps::numeric::real_type<T>::type real(T f){
        return alps::numeric::real(f);
    }

    template<typename _InputIterator, typename _Tp>
    inline _Tp
    accumulate(_InputIterator __first, _InputIterator __last, _Tp __init){
        return std::accumulate(__first, __last, __init);
    }

    template<typename T>
    T sqrt(T arg){
        return std::sqrt(arg);
    }

}
#endif

#endif
