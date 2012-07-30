#ifndef __AMBIENT_NUMERIC_MATRIX_H__
#define __AMBIENT_NUMERIC_MATRIX_H__

#include "ambient/ambient.hpp"
#include "utils/io.hpp"
#include <alps/hdf5.hpp>

#ifndef RVALUE
#define RVALUE
#endif

namespace ambient { namespace numeric {

    template <typename T>
    class weak_matrix_impl;

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
    };

    template <class T>
    class matrix {
    public:
        typedef matrix_impl<T> I;
        typedef typename I::ptr ptr;
        typedef typename I::real_type real_type;
        typedef typename I::size_type size_type; 
        typedef typename I::value_type value_type;
        typedef typename I::scalar_type scalar_type;
        typedef typename I::difference_type difference_type;
        inline void* operator new (size_t);
        inline void operator delete (void*);
        static inline matrix<T> identity_matrix(size_type size);
        explicit inline matrix(const ptr& p);
        explicit inline matrix();
        explicit inline matrix(size_type rows, size_type cols, value_type init_value = value_type()); 
        inline matrix(const matrix& m);
        matrix& operator = (const matrix& rhs); 
#ifdef RVALUE
        inline matrix(matrix&& m); 
        matrix& operator = (matrix&& rhs);
#endif
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
    public:
        ptr impl;
    };

    template <typename T>
    class matrix_impl : public ambient::iteratable<ambient::history> {
    public:
        typedef T value_type;
        typedef size_t size_type;
        typedef ptrdiff_t difference_type;
        typedef typename boost::intrusive_ptr<matrix_impl<T> > ptr;
        typedef typename ambient::future<double> real_type;
        typedef typename ambient::future<T> scalar_type;
        inline void* operator new (size_t);
        inline void operator delete (void*);
        inline matrix_impl();
        inline matrix_impl(size_type rows, size_type cols);
        inline matrix_impl(matrix_impl const& m);
        inline bool empty() const;
        inline size_type num_rows() const; 
        inline size_type num_cols() const;
        inline value_type& get(size_type i, size_type j);
        friend void intrusive_ptr_add_ref(matrix_impl* p){ ++(p->references); } 
        friend void intrusive_ptr_release(matrix_impl* p){ if(--(p->references) == 0) delete p; } 
        operator weak_matrix_impl<T>& ();
    private:
        long references;
    };

    template <typename T>
    class weak_matrix_impl : public matrix_impl<T> { 
        typedef typename boost::intrusive_ptr<matrix_impl<T> > ptr;
    };

} }

#endif
