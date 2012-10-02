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
                    ambient::numeric::kernels::cast_to_vector<T>::spawn(v_ptr, m[k][kk], num_rows, num_cols, num_rows, offset);
                    offset += num_rows;
                }
            }
            ambient::sync();
            return set;
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
                ambient::numeric::kernels::cast_to_vector<T>::spawn(v_ptr, m[k], num_rows, num_cols, num_rows, offset);
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
                    ambient::numeric::kernels::cast_from_vector<T>::spawn(v_ptr, tile, rows, cols, lda, offset);
                    offset += rows;
                }
            }

            ambient::sync();
            return pm;
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
            ambient::numeric::kernels::cast_from_vector<T>::spawn(v_ptr, pm, num_rows, num_cols, lda, offset);
            ambient::sync();
            return pm;
        }
    };

    template <typename T>
    struct binding< alps::numeric::matrix<T>, ambient::numeric::tiles<ambient::numeric::matrix<T> > > {
        static alps::numeric::matrix<T> convert(const ambient::numeric::tiles<ambient::numeric::matrix<T> >& pm){
            size_t num_rows = pm.num_rows();
            size_t num_cols = pm.num_cols();
            size_t offset(0);
            alps::numeric::matrix<T> m(num_rows, num_cols);    
            std::vector<typename alps::numeric::matrix<T>::value_type>* v_ptr = &m.get_values();
            size_t lda = m.stride2();

            for(size_t j = 0; j < pm.nt; ++j){
                size_t offset = j*lda*AMBIENT_IB;
                for(size_t i = 0; i < pm.mt; ++i){
                    const ambient::numeric::matrix<T>& tile = pm.tile(i,j);
                    size_t rows = tile.num_rows();
                    size_t cols = tile.num_cols();
                    ambient::numeric::kernels::cast_to_vector<T>::spawn(v_ptr, tile, rows, cols, lda, offset);
                    offset += rows;
                }
            }

            ambient::sync();
            return m;
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
            ambient::numeric::kernels::cast_to_vector<T>::spawn(v_ptr, pm, num_rows, num_cols, num_rows, offset);
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
            ambient::numeric::diagonal_matrix<T> pm(num_rows, num_rows);    
            const std::vector<typename alps::numeric::diagonal_matrix<T>::value_type>* v_ptr = &m.get_values();
            ambient::numeric::kernels::cast_from_vector<T>::spawn(v_ptr, pm, num_rows, num_cols, num_rows, offset);
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
            alps::numeric::diagonal_matrix<T> m(num_rows, num_rows);    
            std::vector<typename alps::numeric::diagonal_matrix<T>::value_type>* v_ptr = &m.get_values();
            ambient::numeric::kernels::cast_to_vector<T>::spawn(v_ptr, pm, num_rows, num_cols, num_rows, offset);
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
    struct real_type<ambient::future<std::complex<double> > > {
        typedef typename ambient::future<double> type;
    };

    inline double real(double f){
        return f;
    }

    inline double real(const std::complex<double>& f){
        return f.real();
    }

    template<typename T>
    inline const typename real_type<ambient::future<T> >::type& 
    real(const ambient::future<T>& f){
        return ambient::real(f);
    }

    template <typename T,                          // value_type
             template<class AT> class A,           // allocator 
             template<class TT, class AA> class C> // vector<value_type, allocator>
    inline const C<typename real_type<T>::type, A<typename real_type<T>::type> >& 
    real(const C<ambient::future<T>, A<ambient::future<T> > >& f){
        return ambient::real(f);
    }

}
namespace ietl {

    inline double real(const ambient::future<std::complex<double> >& f){
        return ambient::real(f);
    }

}

#else

namespace maquis {

    template <class T> 
    inline typename alps::numeric::real_type<T>::type real(T f){
        return alps::numeric::real(f);
    }

}
#endif

#endif
