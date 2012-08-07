#ifndef __AMBIENT_NUMERIC_MATRIX_HPP__
#define __AMBIENT_NUMERIC_MATRIX_HPP__

#include "ambient/numeric/matrix/matrix.h"
#include "ambient/numeric/matrix/matrix_algorithms.hpp"

namespace ambient { namespace numeric {

    // {{{ transpose_view

    template<class Matrix>
    inline void* transpose_view<Matrix>::operator new (size_t size){
        return boost::singleton_pool<ambient::utils::empty, sizeof(transpose_view<Matrix>)>::malloc(); 
    }

    template<class Matrix>
    inline void transpose_view<Matrix>::operator delete (void* ptr){
        boost::singleton_pool<ambient::utils::empty, sizeof(transpose_view<Matrix>)>::free(ptr); 
    }

    template <class Matrix>
    transpose_view<Matrix>::transpose_view(const Matrix& m)
    : impl(m.impl) 
    { 
    }

    template <class Matrix>
    transpose_view<Matrix>::operator Matrix () const { 
        return transpose(Matrix(this->impl)); 
    }

    template<class Matrix>
    template<class M> 
    size_t transpose_view<Matrix>::rows(const M& m){ 
        return m.num_cols(); 
    } 

    template<class Matrix>
    template<class M> 
    size_t transpose_view<Matrix>::cols(const M& m){ 
        return m.num_rows(); 
    } 

    template<class Matrix>
    const char* transpose_view<Matrix>::code(){
        return "T"; 
    }  

    // }}}

    // {{{ matrix
    #define size_type   typename matrix<T>::size_type
    #define value_type  typename matrix<T>::value_type
    #define scalar_type typename matrix<T>::scalar_type

    template<typename T>
    inline void* matrix<T>::operator new (size_t size){
        return boost::singleton_pool<ambient::utils::empty, sizeof(matrix<T>)>::malloc(); 
    }

    template<typename T>
    inline void* matrix<T>::operator new (size_t size, void* placement){
        return placement; 
    }

    template<typename T>
    inline void matrix<T>::operator delete (void* ptr){
        boost::singleton_pool<ambient::utils::empty, sizeof(matrix<T>)>::free(ptr); 
    }

    template <typename T>
    inline matrix<T> matrix<T>::identity_matrix(size_type size){
        matrix i(size, size);
        fill_identity(i);
        return i;
    }

    template <typename T>
    inline matrix<T>::matrix(const ptr& p, size_t r) 
    : impl(p), ref(r)
    {
    }

    template <typename T>
    inline matrix<T>::matrix(){ 
        this->impl = new I(ambient::dim2(0,0), sizeof(T)); 
    }

    template <typename T>
    inline matrix<T>::matrix(size_type rows, size_type cols, value_type init_value){
        this->impl = new I(ambient::dim2(cols, rows), sizeof(T)); 
        fill_value(*this, init_value);
    }

    template <typename T>
    inline matrix<T>::matrix(const matrix& m){
        this->impl = new I(m.impl->spec.dim, sizeof(T));
        copy(*this, m);
    }
    
    template <typename T>
    matrix<T>& matrix<T>::operator = (const matrix& rhs){
        assert(!rhs.impl->weak());
        matrix c(rhs);
        this->swap(c);
        return *this;
    }

#ifdef RVALUE
    template <typename T>
    inline matrix<T>::matrix(matrix&& m){
        this->impl = m.impl;
    }

    template <typename T>
    matrix<T>& matrix<T>::operator = (matrix&& rhs){
        this->swap(rhs);
        return *this;
    }
#endif

    template<typename T>
    template<class M> 
    size_t matrix<T>::rows(const M& m){ 
        return m.num_rows(); 
    } 

    template<typename T>
    template<class M> 
    size_t matrix<T>::cols(const M& m){ 
        return m.num_cols(); 
    }

    template<typename T>
    inline size_type matrix<T>::num_rows() const { 
        return this->impl->spec.dim.y; 
    }

    template<typename T>
    inline size_type matrix<T>::num_cols() const {
        return this->impl->spec.dim.x; 
    }

    template<typename T>
    inline scalar_type matrix<T>::trace() const { 
        return trace(*this);           
    }

    template<typename T>
    inline void matrix<T>::transpose(){ 
        transpose_inplace(*this);      
    }

    template<typename T>
    inline void matrix<T>::conj(){ 
        conj_inplace(*this);           
    }

    template<typename T>
    inline bool matrix<T>::empty() const { 
        return (this->impl->spec.dim == 0);    
    }

    template<typename T>
    inline void matrix<T>::swap(matrix& r){ 
        this->impl.swap(r.impl);       
    }

    template<typename T>
    inline void matrix<T>::resize(size_type rows, size_type cols){
        ambient::numeric::resize(*this, rows, cols);
    }

    template<typename T>
    inline void matrix<T>::remove_rows(size_type i, size_type k){
        remove_rows(*this, i, k);
    }

    template<typename T>
    inline void matrix<T>::remove_cols(size_type j, size_type k){
        remove_cols(*this, j, k); 
    }

    template<typename T>
    inline matrix<T>& matrix<T>::operator += (const matrix& rhs){
        add_inplace(*this, rhs);
        return *this;
    }

    template<typename T>
    inline matrix<T>& matrix<T>::operator -= (const matrix& rhs){
        sub_inplace(*this, rhs);
        return *this;
    }

    template<typename T>
    template <typename T2> 
    inline matrix<T>& matrix<T>::operator *= (const T2& t){
        mul_inplace(*this, t);
        return *this;
    }

    template<typename T>
    template <typename T2> 
    inline matrix<T>& matrix<T>::operator /= (const T2& t){
        div_inplace(*this, t);
        return *this;
    }

    template<typename T>
    inline value_type& matrix<T>::operator() (size_type i, size_type j){
        ambient::playout(); return ((value_type*)ambient::controller.ufetch(*this->impl->current))[ j*this->impl->spec.dim.y + i ];
    }

    template<typename T>
    inline const value_type& matrix<T>::operator() (size_type i, size_type j) const {
        ambient::playout(); return ((value_type*)ambient::controller.ufetch(*this->impl->current))[ j*this->impl->spec.dim.y + i ];
    }

    template<typename T>
    const char* matrix<T>::code(){ 
        return "N"; 
    }  
    #undef size_type
    #undef value_type
    #undef scalar_type
    // }}}

} }

#endif
