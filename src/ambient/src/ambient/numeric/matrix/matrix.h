#ifndef __AMBIENT_NUMERIC_MATRIX_H__
#define __AMBIENT_NUMERIC_MATRIX_H__

#include "ambient/ambient.hpp"
#include "types/utils/traits.hpp"
#include "utils/io.hpp"

#ifdef HAVE_ALPS_HDF5
#include <alps/hdf5.hpp>
#endif

namespace ambient { namespace numeric {

    template<class T>
    class diagonal_matrix;

    template <class T>
    class matrix {
    public:
        typedef matrix_impl<T> I;
        typedef typename I::ptr ptr;
        typedef typename I::value_type value_type;
        typedef typename I::scalar_type scalar_type;
        typedef typename I::size_type size_type; 
        typedef typename I::difference_type difference_type;
        // {{{ matrix_impl forwarding

        static inline matrix<T> identity_matrix(size_type size){
            return matrix<T>(size);
        }

        inline matrix(size_type size){ // identity
            this->impl = new I(size, size); 
            this->impl->fill_identity();
        }

        inline matrix(){ // shouldn't be called
            this->impl = new I(); 
        }

        inline matrix(size_type rows, size_type cols, value_type init_value = value_type()){
            this->impl = new I(rows, cols); 
            this->impl->fill_value(init_value);
        }

        inline matrix(const matrix& m){
            this->impl = new I(*m.impl);
            this->impl->copy(*m.impl);
        }

        matrix& operator = (matrix rhs){
            assert(!rhs.impl->weak());
            this->swap(rhs);
            return *this;
        }

    public:

        inline void swap(matrix& r){
            this->impl.swap(r.impl);
        }

        friend void swap(matrix& x, matrix& y){
            x.swap(y);
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

        inline void conj(){
            this->impl->conj();
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
            assert(this->num_rows() != 0 && this->num_cols() != 0);
            if(this->num_rows() == rows && this->num_cols() == cols) return;
            matrix resized(rows, cols);
            if(!this->impl->weak())
                this->impl->resize(*resized.impl, rows, cols);
            this->swap(resized);
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

        inline matrix& operator += (const matrix& rhs){
            this->impl->add(*rhs.impl);
            return *this;
        }

        inline matrix& operator -= (const matrix& rhs){
            this->impl->sub(*rhs.impl);
            return *this;
        }

        inline matrix& operator *= (const matrix& rhs){
            this->impl->mul(*rhs.impl);
            return *this;
        }

        inline matrix& operator *= (const diagonal_matrix<T>& rhs){
            this->impl->mul(rhs);
            return *this;
        }

        template <typename T2> inline matrix& operator *= (const T2& t){
            this->impl->mul(t);
            return *this;
        }

        template <typename T2> inline matrix& operator /= (const T2& t){
            this->impl->div(t);
            return *this;
        }

#ifdef HAVE_ALPS_HDF5
        inline void load(alps::hdf5::archive & ar) { /*ambient::cerr << "I don't do much." << std::endl;*/ }
        inline void save(alps::hdf5::archive & ar) const { /*ambient::cerr << "I don't do much either." << std::endl;*/ }
#endif
        // }}}
    public:
        ptr impl;
    };

    template <typename T>
    class matrix_impl :
    public ambient::iteratable<ambient::history>
    {
    public:
        typedef T         value_type;      // The type T of the elements of the matrix
        typedef size_t    size_type;       // Unsigned integer type that represents the dimensions of the matrix
        typedef ptrdiff_t difference_type; // Signed integer type to represent the distance of two elements in the memory
        typedef typename boost::intrusive_ptr<matrix_impl<T> > ptr;
        typedef typename ambient::future<T> scalar_type;

        inline matrix_impl();             // please avoid implicit conversions
        inline matrix_impl(size_type rows, size_type cols);
        inline matrix_impl(matrix_impl const& m);
        inline value_type& get(size_type i, size_type j);
        inline scalar_type trace() const;
        inline void fill_identity();
        inline void fill_random();
        inline void fill_value(value_type v);
        inline void conj();
        inline void transpose();
        inline bool empty() const;
        inline size_type num_rows() const;
        inline size_type num_cols() const;
        inline bool atomic() const;
        inline void resize(matrix_impl& r, size_type rows, size_type cols);
        inline void remove_rows(size_type i, size_type k);
        inline void remove_cols(size_type j, size_type k);
        inline void add(const matrix_impl& rhs); 
        inline void sub(const matrix_impl& rhs);
        inline void mul(const matrix_impl& rhs);
        inline void mul(const diagonal_matrix<T>& rhs);
        template <typename T2> inline void mul(const T2& t);
        template <typename T2> inline void div(const T2& t);
        inline void copy(const matrix_impl& m);
        friend void intrusive_ptr_add_ref(matrix_impl* p){ ++(p->references); }
        friend void intrusive_ptr_release(matrix_impl* p){ if(--(p->references) == 0) delete p; }
    private:
        long references;
    };

} } // namespace ambient::numeric

#endif
