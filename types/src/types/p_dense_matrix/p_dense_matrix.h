#ifndef __MAQUIS_TYPES_P_DENSE_MATRIX_HPP__
#define __MAQUIS_TYPES_P_DENSE_MATRIX_HPP__

#include "ambient/ambient.hpp"
#include "types/utils/traits.hpp"
#include "utils/io.hpp"

#ifdef HAVE_ALPS_HDF5
#include <alps/hdf5.hpp>
#endif

namespace maquis { namespace types {

    template<class T>
    class p_diagonal_matrix;

    template <class T>
    class p_dense_matrix {
    public:
        typedef p_dense_matrix_impl<T> I;
        typedef typename I::ptr ptr;
        typedef typename I::value_type value_type;
        typedef typename I::scalar_type scalar_type;
        typedef typename I::size_type size_type; 
        typedef typename I::difference_type difference_type;
        // {{{ p_dense_matrix_impl forwarding

        inline p_dense_matrix(){
            this->impl = new I(); 
        }

        inline p_dense_matrix(size_type rows, size_type cols, value_type init_value = value_type()){
            this->impl = new I(rows, cols); 
            this->impl->fill_value(init_value);
        }

        inline p_dense_matrix(p_dense_matrix const& m){
            this->impl = new I(*m.impl);
        }

        static inline p_dense_matrix<T> identity_matrix(size_type size){
            p_dense_matrix<T> m(size, size);
            m.fill_identity();
            return m;
        }

        inline value_type& get(size_type i, size_type j) const {
            return this->impl->get(i,j);
        }

        inline scalar_type trace() const {
            return this->impl->trace();
        }

        inline void fill_identity(){
            this->impl->fill_identity();
        }

        inline void fill_value(value_type v){
            this->impl->fill_value(v);
        }

        inline void fill_random(){
            this->impl->fill_random();
        }

        inline void conjugate(){
            this->impl->conjugate();
        }

        inline void transpose(){
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

        inline bool atomic() const {
            return this->impl->atomic();
        }

        inline void resize(size_type rows, size_type cols){
            if(this->num_rows() != rows || this->num_cols() != cols){
                assert(this->num_rows() != 0 && this->num_cols() != 0);
                p_dense_matrix resized(rows, cols);
                this->impl->resize(*resized.impl, rows, cols); 
                this->impl.swap(resized.impl);
            }
        }

        inline void remove_rows(size_type i, size_type k = 1){
            this->impl->remove_rows(i, k);
            this->resize(this->num_rows()-k, this->num_cols());
        }

        inline void remove_cols(size_type j, size_type k = 1){
            this->impl->remove_cols(j, k); 
            this->resize(this->num_rows(), this->num_cols()-k);
        }

        inline value_type& operator() (size_type i, size_type j){
            return this->impl->get(i,j);
        }

        inline const value_type& operator() (size_type i, size_type j) const {
            return this->impl->get(i,j);
        }

        inline p_dense_matrix& operator = (const p_dense_matrix& rhs){
            this->resize(rhs.num_rows(), rhs.num_cols());
            if(!rhs.impl->pt_clean()) this->impl->cpy(*rhs.impl);
            return *this;
        }

        inline p_dense_matrix& operator += (const p_dense_matrix& rhs){
            this->impl->add(*rhs.impl);
            return *this;
        }

        inline p_dense_matrix& operator -= (const p_dense_matrix& rhs){
            this->impl->sub(*rhs.impl);
            return *this;
        }

        inline p_dense_matrix& operator *= (const p_dense_matrix& rhs){
            this->impl->mul(*rhs.impl);
            return *this;
        }

        inline p_dense_matrix& operator *= (const p_diagonal_matrix<T>& rhs){
            this->impl->mul(rhs);
            return *this;
        }

        template <typename T2> inline p_dense_matrix& operator *= (const T2& t){
            this->impl->mul(t);
            return *this;
        }

        template <typename T2> inline p_dense_matrix& operator /= (const T2& t){
            this->impl->mul(((value_type)1)/t);
            return *this;
        }

#ifdef HAVE_ALPS_HDF5
        inline void load(alps::hdf5::archive & ar) { /*maquis::cerr << "I don't do much." << std::endl;*/ }
        inline void save(alps::hdf5::archive & ar) const { /*maquis::cerr << "I don't do much either." << std::endl;*/ }
#endif
        // }}}
    public:
        ptr impl;
    };

    template <typename T>
    class p_dense_matrix_impl :
    public ambient::parallel< p_dense_matrix_impl<T> >
    {
    public:
        typedef T         value_type;      // The type T of the elements of the matrix
        typedef size_t    size_type;       // Unsigned integer type that represents the dimensions of the matrix
        typedef ptrdiff_t difference_type; // Signed integer type to represent the distance of two elements in the memory
        typedef typename ambient::parallel< p_dense_matrix_impl<T> >::ptr ptr;
        typedef typename ambient::future<T> scalar_type;

        inline ~p_dense_matrix_impl();
        inline p_dense_matrix_impl();             // please avoid implicit conversions
        inline p_dense_matrix_impl(size_type rows, size_type cols);
        inline p_dense_matrix_impl(p_dense_matrix_impl const& m);
        inline value_type& get(size_type i, size_type j);
        inline scalar_type trace() const;
        inline void fill_identity();
        inline void fill_random();
        inline void fill_value(value_type v);
        inline void conjugate();
        inline void transpose();
        inline bool empty() const;
        inline size_type num_rows() const;
        inline size_type num_cols() const;
        inline bool atomic() const;
        inline void resize(p_dense_matrix_impl& r, size_type rows, size_type cols);
        inline void remove_rows(size_type i, size_type k);
        inline void remove_cols(size_type j, size_type k);
        inline void cpy(const p_dense_matrix_impl& rhs);
        inline void add(const p_dense_matrix_impl& rhs); 
        inline void sub(const p_dense_matrix_impl& rhs);
        inline void mul(const p_dense_matrix_impl& rhs);
        inline void mul(const p_diagonal_matrix<T>& rhs);
        template <typename T2> inline void mul(const T2& t);
    private:
        size_type rows;
        size_type cols;
    };

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
        //typedef p_diagonal_matrix<T> type;
        typedef std::vector<T> type;
    };
    
    template<typename T>
    struct associated_real_vector<p_dense_matrix<T> > {
        //typedef p_diagonal_matrix<typename detail::real_type<T>::type> type;
        typedef std::vector<typename detail::real_type<T>::type> type;
    };
    // }}}

} } // namespace maquis::types
#include "types/p_dense_matrix/p_dense_matrix.hpp"
#endif
