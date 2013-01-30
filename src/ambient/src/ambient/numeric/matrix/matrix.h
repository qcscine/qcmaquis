#ifndef AMBIENT_NUMERIC_MATRIX_H
#define AMBIENT_NUMERIC_MATRIX_H

#include "ambient/ambient.hpp"
#include "utils/io.hpp"
#include <alps/hdf5.hpp>

#ifndef RVALUE
#define RVALUE
#endif

namespace ambient { namespace numeric {

    template <class T>
    class weak_view;

    template <typename T>
    class matrix {
    public:
        typedef ambient::history I;
        typedef T value_type;
        typedef size_t size_type;
        typedef ptrdiff_t difference_type;
        typedef typename ambient::history* ptr;
        typedef typename ambient::future<double> real_type;
        typedef typename ambient::future<T> scalar_type;

        void* operator new (size_t);
        void* operator new (size_t, void*);
        void operator delete (void*);
        void operator delete (void*, void*){ } // doesn't throw
        static matrix<T> identity_matrix(size_type size);
        explicit matrix(const ptr& p, size_t r);
        explicit matrix();
        explicit matrix(size_type rows, size_type cols, value_type init_value = value_type()); 
        matrix(const matrix& a);
        matrix& operator = (const matrix& rhs); 
#if 0
        matrix(matrix&& a); 
        matrix& operator = (matrix&& rhs);
#endif
       ~matrix();
    public:
        template<class M> static size_t inc (const M& a); 
        template<class M> static size_t rows(const M& a); 
        template<class M> static size_t cols(const M& a);
        size_type lda() const;
        size_type num_rows() const;
        size_type num_cols() const;
        scalar_type trace() const;
        void transpose();
        void conj();
        bool empty() const;          
        void swap(matrix& r);
        template<typename TT> 
        friend void swap(matrix& x, matrix& y);
        void resize(size_type m, size_type n); 
        void remove_rows(size_type i, size_type k = 1);
        void remove_cols(size_type j, size_type k = 1);
        matrix& locate(size_type i, size_type j);
        const matrix& locate(size_type i, size_type j) const;
        size_t addr(size_type i, size_type j) const;
        matrix& operator += (const matrix& rhs);
        matrix& operator -= (const matrix& rhs);
        template <typename T2> matrix& operator *= (const T2& t);
        template <typename T2> matrix& operator /= (const T2& t);
        value_type& operator() (size_type i, size_type j);
        const value_type& operator() (size_type i, size_type j) const;
        void load(alps::hdf5::archive & ar){};
        void save(alps::hdf5::archive & ar)const{};
        static const char* code();
        operator weak_view<T>& (){ return *(weak_view<T>*)this; }
    public:
        ptr impl;
        size_t ref;
    };

    template <class Matrix>
    class subset_view {
    public:
        typedef typename Matrix::real_type real_type;
        typedef typename Matrix::size_type size_type; 
        typedef typename Matrix::value_type value_type;
        typedef typename Matrix::scalar_type scalar_type;
        typedef typename Matrix::difference_type difference_type;
        typedef typename Matrix::ptr ptr;
        subset_view(const Matrix& a) : impl(a.impl), m(&a) {}
        size_t num_rows(){ return m->num_rows(); };
        size_t num_cols(){ return m->num_cols(); };
        template<class M> static size_t rows(const M& a); 
        template<class M> static size_t cols(const M& a);
        static const char* code();
        operator Matrix& () const { return *(Matrix*)m; }
        ptr impl;
        const Matrix* m;
    };

    template <class T>
    class weak_view : public matrix<T> {
    public:
        weak_view(const typename matrix<T>::ptr& p, size_t r) : matrix<T>(p, r) {}
    };

    template <class Matrix>
    class transpose_view {
    public:
        typedef typename Matrix::real_type real_type;
        typedef typename Matrix::size_type size_type; 
        typedef typename Matrix::value_type value_type;
        typedef typename Matrix::scalar_type scalar_type;
        typedef typename Matrix::difference_type difference_type;
        void* operator new (size_t);
        void operator delete (void*);
        explicit transpose_view(const Matrix& a);
        transpose_view& locate(size_type i, size_type j);
        const transpose_view& locate(size_type i, size_type j) const;
        size_t addr(size_type i, size_type j) const;
        size_t lda() const;
        operator Matrix () const;
        template<class M> static size_t inc (const M& a); 
        template<class M> static size_t rows(const M& a); 
        template<class M> static size_t cols(const M& a);
        static const char* code();
        typename Matrix::ptr impl;
        size_t ref;
    };

} }

#endif
