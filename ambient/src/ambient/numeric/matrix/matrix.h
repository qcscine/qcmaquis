#ifndef __AMBIENT_NUMERIC_MATRIX_H__
#define __AMBIENT_NUMERIC_MATRIX_H__

#include "ambient/ambient.hpp"
#include "utils/io.hpp"
#include <alps/hdf5.hpp>

#ifndef RVALUE
#define RVALUE
#endif

namespace ambient { namespace numeric {

    template <class T>
    class weak_view;

    template <class Matrix>
    class transpose_view {
    public:
        typedef typename Matrix::real_type real_type;
        typedef typename Matrix::size_type size_type; 
        typedef typename Matrix::value_type value_type;
        typedef typename Matrix::scalar_type scalar_type;
        inline void* operator new (size_t);
        inline void operator delete (void*);
        explicit transpose_view(const Matrix& m);
        operator Matrix () const;
        template<class M> static size_t rows(const M& m); 
        template<class M> static size_t cols(const M& m);
        static const char* code();
        typename Matrix::ptr impl;
        size_t ref;
    };

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

        inline void* operator new (size_t);
        inline void* operator new (size_t, void*);
        inline void operator delete (void*);
        static inline matrix<T> identity_matrix(size_type size);
        explicit inline matrix(const ptr& p, size_t r);
        explicit inline matrix();
        explicit inline matrix(size_type rows, size_type cols, value_type init_value = value_type()); 
        inline matrix(const matrix& m);
        matrix& operator = (const matrix& rhs); 
#ifdef RVALUE
        inline matrix(matrix&& m); 
        matrix& operator = (matrix&& rhs);
#endif
        inline ~matrix();
    public:
        template<class M> static size_t rows(const M& m); 
        template<class M> static size_t cols(const M& m);
        inline size_type num_rows() const;
        inline size_type num_cols() const;
        inline scalar_type trace() const;
        inline void transpose();
        inline void conj();
        inline bool empty() const;          
        inline void swap(matrix& r);
        friend void swap(matrix& x, matrix& y);
        inline void resize(size_type rows, size_type cols); 
        inline void remove_rows(size_type i, size_type k = 1);
        inline void remove_cols(size_type j, size_type k = 1);
        inline matrix& operator += (const matrix& rhs);
        inline matrix& operator -= (const matrix& rhs);
        template <typename T2> inline matrix& operator *= (const T2& t);
        template <typename T2> inline matrix& operator /= (const T2& t);
        inline value_type& operator() (size_type i, size_type j);
        inline const value_type& operator() (size_type i, size_type j) const;
        inline void load(alps::hdf5::archive & ar){};
        inline void save(alps::hdf5::archive & ar)const{};
        static const char* code();
        operator weak_view<T>& (){ return *(weak_view<T>*)this; }
    public:
        ptr impl;
        size_t ref;
    };

    template <class T>
    class weak_view : public matrix<T> {
        public:
        operator matrix<T>& (){ return *(matrix<T>*)this; }
        weak_view(const typename matrix<T>::ptr& p, size_t r) : matrix<T>(p, r) {}
    };

} }

#endif
