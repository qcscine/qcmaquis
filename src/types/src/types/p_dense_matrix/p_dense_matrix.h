#ifndef __MAQUIS_TYPES_P_DENSE_MATRIX_HPP__
#define __MAQUIS_TYPES_P_DENSE_MATRIX_HPP__

#include "ambient/ambient.hpp"
#include "types/utils/traits.hpp"

#ifdef HDF5
#include <alps/hdf5.hpp>
#endif

#include <ostream>

#include "types/p_dense_matrix/kernels/i_kernels.hpp"

namespace maquis { namespace types {

    template <class T>
    class p_dense_matrix {
    public:
        typedef p_dense_matrix_impl<T> I;
        typedef typename I::ptr ptr;
        typedef typename I::value_type value_type;
        typedef typename I::size_type size_type; 
        typedef typename I::difference_type difference_type;
        // {{{ p_dense_matrix_impl forwarding
        p_dense_matrix(){ 
            this->impl = new I(); 
        }

        p_dense_matrix(size_type rows, size_type cols = 0, value_type init_value = value_type()){
            this->impl = new I(rows, cols, init_value); 
        }

        p_dense_matrix(p_dense_matrix const& m){
            this->impl = new I(*m.impl);
        }

        static p_dense_matrix<T> identity_matrix(size_type size){
            p_dense_matrix m(size, size);
            m.impl->pt_set_init(ambient::identity_i<T>);
            return m;
        }

        void generate(){
            this->impl->pt_set_init(ambient::random_i<T>);
        }

        template<typename O>
        void set_init(void(*fp)(O&)){
            this->impl->pt_set_init(fp);
        }

        value_type& get(size_type i, size_type j) const {
            return this->impl->get(i,j);
        }

        value_type trace(){
            return this->impl->trace();
        }

        void conjugate(){
            this->impl->conjugate();
        }

        void transpose(){
            this->impl->transpose();
        }

        inline bool empty() const {
            return this->impl->empty();
        }

        inline size_type num_rows() const {
            return this->impl->num_rows();
        }

        inline size_type num_cols() const {
            return this->impl->num_cols();
        }

        void resize(size_type rows, size_type cols){
            this->impl->resize(rows, cols); 
        }

        void remove_rows(size_type i, size_type k){
            this->impl->remove_rows(i, k);
        }

        void remove_cols(size_type j, size_type k){
            this->impl->remove_cols(j, k); 
        }

        void clear(){
            this->impl->clear();
        }

        inline value_type& operator() (size_type i, size_type j){
            return this->impl->get(i,j);
        }

        inline const value_type& operator() (size_type i, size_type j) const {
            return this->impl->get(i,j);
        }

        p_dense_matrix& operator = (const p_dense_matrix& rhs){
            this->impl->cpy(*rhs.impl);
            return *this;
        }

        p_dense_matrix& operator += (const p_dense_matrix& rhs){
            this->impl->add(*rhs.impl);
            return *this;
        }

        p_dense_matrix& operator -= (const p_dense_matrix& rhs){
            this->impl->sub(*rhs.impl);
            return *this;
        }

        p_dense_matrix& operator *= (const p_dense_matrix& rhs){
            this->impl->mul(*rhs.impl);
            return *this;
        }

        template <typename T2> p_dense_matrix& operator *= (const T2& t){
            this->impl->mul(t);
            return *this;
        }

        template <typename T2> p_dense_matrix& operator /= (const T2& t){
            this->impl->mul(1/t);
            return *this;
        }

#ifdef HAVE_ALPS_HDF5
        void load(alps::hdf5::archive & ar) { maquis::cerr << "I don't do much." << std::endl; }
        void save(alps::hdf5::archive & ar) const { maquis::cerr << "I don't do much either." << std::endl; }
#endif
        // }}}
    public:
        ptr impl;
    };

    // {{{ p_dense_matrix operators (free functions)
    #define size_type typename ::maquis::types::p_dense_matrix_impl<T>::size_type
    template <typename T>
    const p_dense_matrix<T>& operator + (p_dense_matrix<T> lhs, const p_dense_matrix<T>& rhs){
        return (lhs += rhs); 
    }

    template <typename T>
    const p_dense_matrix<T>& operator - (p_dense_matrix<T> lhs, const p_dense_matrix<T>& rhs){ 
        return (lhs -= rhs); 
    }

    template<typename T>
    const p_dense_matrix<T>& operator * (p_dense_matrix<T> lhs, const p_dense_matrix<T>& rhs){ 
        return (lhs *= rhs); 
    }

    template<typename T, typename T2>
    const p_dense_matrix<T>& operator * (p_dense_matrix<T> lhs, const T2& rhs){ 
        return (lhs *= rhs); 
    }

    template<typename T, typename T2>
    const p_dense_matrix<T>& operator * (const T2& lhs, p_dense_matrix<T> rhs){ 
        return (rhs *= lhs); 
    }

    template <typename T>
    std::ostream& operator << (std::ostream& o, p_dense_matrix<T> const& m){
        if(ambient::outlet())
        for(size_type i=0; i< m.num_rows(); ++i){
            for(size_type j=0; j < m.num_cols(); ++j)
                printf("%.4f	", m(i,j));
            printf("\n");
        }
        return o;
    }

    template<typename T>
    size_type num_rows(const p_dense_matrix<T>& m){
        return m.num_rows();
    }

    template<typename T>
    size_type num_cols(const p_dense_matrix<T>& m){
        return m.num_cols();
    }

    template<typename T>
    void resize(p_dense_matrix<T>& m, size_t rows, size_t cols){
        m.resize(rows, cols);
    }

    template<typename T>
    T trace(p_dense_matrix<T>& m){
        return m.trace();
    }

    template<typename T>
    p_dense_matrix<T> transpose(p_dense_matrix<T> m){
        m.transpose();
        return m;
    }
    #undef size_type
    // }}}

    // p_dense_matrix realization (derives from parallel_t) //
    template <typename T>
    class p_dense_matrix_impl :
    public ambient::parallel_t< p_dense_matrix_impl<T> >
    {
    private:
        p_dense_matrix_impl();                            // avoiding implicit conversions
    public:
        typedef T         value_type;                     // The type T of the elements of the matrix
        typedef size_t    size_type;                      // Unsigned integer type that represents the dimensions of the matrix
        typedef ptrdiff_t difference_type;                // Signed integer type to represent the distance of two elements in the memory
        typedef typename ambient::parallel_t< p_dense_matrix_impl<T> >::ptr ptr;

       ~p_dense_matrix_impl();
        p_dense_matrix_impl(size_type rows, size_type cols, T init_value);
        p_dense_matrix_impl(p_dense_matrix_impl const& m);
        value_type& get(size_type i, size_type j);
        value_type trace();
  
        void conjugate();
        void transpose();
        inline bool empty() const;
        inline size_type num_rows() const;
        inline size_type num_cols() const;
        void resize(size_type rows, size_type cols);
        void remove_rows(size_type i, size_type k);
        void remove_cols(size_type j, size_type k);
        void clear();

        void cpy(const p_dense_matrix_impl& rhs);
        void add(const p_dense_matrix_impl& rhs); 
        void sub(const p_dense_matrix_impl& rhs);
        void mul(const p_dense_matrix_impl& rhs);
        template <typename T2> void mul(const T2& t);
    private:
        value_type init_value;
        size_type rows;
        size_type cols;
    };

    template<class T>
    class p_diagonal_matrix;

    // {{{ matrix-specific associated types
    template<typename T>
    struct associated_diagonal_matrix< p_dense_matrix<T> > {
        typedef p_diagonal_matrix<T> type;
    };

    template<typename T>
    struct associated_real_diagonal_matrix< p_dense_matrix<T> > {
        typedef p_diagonal_matrix<typename detail::real_type<T>::type> type;
    };

    template<typename T>
    struct associated_vector<p_dense_matrix<T> > {
        typedef p_diagonal_matrix<T> type;
    };
    
    template<typename T>
    struct associated_real_vector<p_dense_matrix<T> > {
        typedef p_diagonal_matrix<typename detail::real_type<T>::type> type;
    };
    // }}}
} } // namespace maquis::types
#include "types/p_dense_matrix/p_dense_matrix_impl.hpp"
#endif
