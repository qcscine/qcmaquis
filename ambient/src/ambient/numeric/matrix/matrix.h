/*
 * Ambient, License - Version 1.0 - May 3rd, 2012
 *
 * Permission is hereby granted, free of charge, to any person or organization
 * obtaining a copy of the software and accompanying documentation covered by
 * this license (the "Software") to use, reproduce, display, distribute,
 * execute, and transmit the Software, and to prepare derivative works of the
 * Software, and to permit third-parties to whom the Software is furnished to
 * do so, all subject to the following:
 *
 * The copyright notices in the Software and this entire statement, including
 * the above license grant, this restriction and the following disclaimer,
 * must be included in all copies of the Software, in whole or in part, and
 * all derivative works of the Software, unless such copies or derivative
 * works are solely in the form of machine-executable object code generated by
 * a source language processor.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

#ifndef AMBIENT_NUMERIC_MATRIX_H
#define AMBIENT_NUMERIC_MATRIX_H

#include "ambient/ambient.hpp"
#include "ambient/numeric/future.hpp"

namespace ambient { namespace numeric {

    template <typename T, class Allocator = ambient::default_allocator<T> >
    class matrix : public ambient::memory::use_fixed_new<matrix<T,Allocator> > {
    public:
        typedef T value_type;
        typedef size_t size_type;
        typedef Allocator allocator_type;
        typedef size_t difference_type;
        typedef typename ambient::numeric::future<double> real_type;
        typedef typename ambient::numeric::future<T> scalar_type;

        explicit matrix();
        explicit matrix(size_type rows, size_type cols, value_type init_value = value_type()); 
        matrix(const matrix& a);
        matrix& operator = (const matrix& rhs);
        template<class OtherAllocator>
        matrix& operator = (const matrix<T,OtherAllocator>& rhs); 
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
        void swap(matrix& r);
        template<typename TT> 
        friend void swap(matrix& x, matrix& y);
        void resize(size_type m, size_type n); 
        matrix& locate(size_type i, size_type j);
        const matrix& locate(size_type i, size_type j) const;
        size_t addr(size_type i, size_type j) const;
        matrix& operator += (const matrix& rhs);
        matrix& operator -= (const matrix& rhs);
        template <typename T2> matrix& operator *= (const T2& t);
        template <typename T2> matrix& operator /= (const T2& t);
        value_type& operator() (size_type i, size_type j);
        const value_type& operator() (size_type i, size_type j) const;
        static const char* code();
    public:
        ambient_version( T data[ AMBIENT_VAR_LENGTH ]; );
    };

    template <class Matrix>
    class subset_view {
    public:
        typedef typename Matrix::real_type real_type;
        typedef typename Matrix::size_type size_type; 
        typedef typename Matrix::value_type value_type;
        typedef typename Matrix::scalar_type scalar_type;
        typedef typename Matrix::difference_type difference_type;
        typedef typename Matrix::allocator_type allocator_type;
        subset_view(const Matrix& a) : versioned(a.versioned), m(&a) {}
        size_t num_rows(){ return m->num_rows(); };
        size_t num_cols(){ return m->num_cols(); };
        template<class M> static size_t rows(const M& a); 
        template<class M> static size_t cols(const M& a);
        static const char* code();
        operator Matrix& () const { return *(Matrix*)m; }
        ambient_non_destroyable ambient_version( value_type data[ AMBIENT_VAR_LENGTH ]; );
        const Matrix* m;
    };

    template <class Matrix>
    class transpose_view : public ambient::memory::use_fixed_new<transpose_view<Matrix> > {
    public:
        typedef typename Matrix::real_type real_type;
        typedef typename Matrix::size_type size_type; 
        typedef typename Matrix::value_type value_type;
        typedef typename Matrix::scalar_type scalar_type;
        typedef typename Matrix::difference_type difference_type;
        typedef typename Matrix::allocator_type allocator_type;
        explicit transpose_view(const Matrix& a) : versioned(a.versioned) {}
        transpose_view& locate(size_type i, size_type j);
        const transpose_view& locate(size_type i, size_type j) const;
        size_t addr(size_type i, size_type j) const;
        size_t lda() const;
        operator Matrix () const;
        template<class M> static size_t inc (const M& a); 
        template<class M> static size_t rows(const M& a); 
        template<class M> static size_t cols(const M& a);
        static const char* code();
        ambient_non_destroyable ambient_version( value_type data[ AMBIENT_VAR_LENGTH ]; );
    };

} }

#endif
