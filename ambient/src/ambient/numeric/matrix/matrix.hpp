#ifndef AMBIENT_NUMERIC_MATRIX_HPP
#define AMBIENT_NUMERIC_MATRIX_HPP

#include "ambient/numeric/matrix/matrix.h"
#include "ambient/numeric/matrix/matrix_algorithms.hpp"

namespace ambient { namespace numeric {

    // {{{ transpose_view

    template<class Matrix>
    inline void* transpose_view<Matrix>::operator new (size_t size){
        return ambient::pool.malloc<transpose_view<Matrix> >(); 
    }

    template<class Matrix>
    inline void transpose_view<Matrix>::operator delete (void* ptr){
        ambient::pool.free<transpose_view<Matrix> >(ptr); 
    }

    template <class Matrix>
    transpose_view<Matrix>::transpose_view(const Matrix& a)
    : core(a.core) 
    { 
    }

    template <class Matrix>
    inline transpose_view<Matrix>& transpose_view<Matrix>::locate(size_type i, size_type j){
        return *this;
    }

    template <class Matrix>
    inline const transpose_view<Matrix>& transpose_view<Matrix>::locate(size_type i, size_type j) const {
        return *this;
    }

    template<class Matrix>
    inline size_t transpose_view<Matrix>::lda() const { 
        return this->core->spec.dim.y; 
    }

    template <class Matrix>
    inline size_t transpose_view<Matrix>::addr(size_type i, size_type j) const {
        return (j + i*lda()); 
    }

    template <class Matrix>
    transpose_view<Matrix>::operator Matrix () const {
        Matrix t(Matrix(this->core,0));
        transpose_inplace(t); 
        return t;
    }

    template<class Matrix>
    template<class M> 
    size_t transpose_view<Matrix>::inc(const M& a){ 
        return a.num_rows(); 
    } 

    template<class Matrix>
    template<class M> 
    size_t transpose_view<Matrix>::rows(const M& a){ 
        return a.num_cols(); 
    } 

    template<class Matrix>
    template<class M> 
    size_t transpose_view<Matrix>::cols(const M& a){ 
        return a.num_rows(); 
    } 

    template<class Matrix>
    const char* transpose_view<Matrix>::code(){
        return "T"; 
    }  

    // }}}
    // {{{ subset_view

    template<class Matrix>
    template<class M> 
    size_t subset_view<Matrix>::rows(const M& a){ 
        return ((Matrix&)a).num_rows(); 
    } 

    template<class Matrix>
    template<class M> 
    size_t subset_view<Matrix>::cols(const M& a){ 
        return ((Matrix&)a).num_cols(); 
    } 

    template<class Matrix>
    const char* subset_view<Matrix>::code(){
        return "N"; 
    }  
    
    // }}}
    // {{{ matrix
    #define size_type   typename matrix<T>::size_type
    #define value_type  typename matrix<T>::value_type
    #define scalar_type typename matrix<T>::scalar_type

    template<typename T>
    inline void* matrix<T>::operator new (size_t size){
        return ambient::pool.malloc<matrix<T> >();
    }

    template<typename T>
    inline void* matrix<T>::operator new (size_t size, void* placement){
        return placement; 
    }

    template<typename T>
    inline void matrix<T>::operator delete (void* ptr){
        ambient::pool.free<matrix<T> >(ptr);
    }

    template <typename T>
    inline matrix<T>::~matrix(){
        if(this->core->weak()) delete this->core;
        else ambient::destroy(this->core);
    }

    template <typename T>
    inline matrix<T>::matrix(const ptr& p, size_t r) 
    : core(p), ref(r)
    {
    }

    template <typename T>
    inline matrix<T>::matrix(){ 
        this->core = new I(ambient::dim2(0,0), sizeof(T)); 
    }

    template <typename T>
    inline matrix<T>::matrix(size_type rows, size_type cols, value_type init_value){
        this->core = new I(ambient::dim2(cols, rows), sizeof(T));
        fill_value(*this, init_value);
    }

    template <typename T>
    inline matrix<T>::matrix(const matrix& a){
        this->core = new I(a.core->spec.dim, sizeof(T));
        ambient::fuse(a.core, this->core);
    }
    
    template <typename T>
    matrix<T>& matrix<T>::operator = (const matrix& rhs){
        //assert(!rhs.core->weak()); // can be optimized if weak
        matrix c(rhs);
        this->swap(c);
        return *this;
    }

#if 0
    template <typename T>
    inline matrix<T>::matrix(matrix&& a){
        printf("ERROR: NOT TESTED (RVALUE COPY)\n");
        this->core = a.core; // need to clear a.core
    }

    template <typename T>
    matrix<T>& matrix<T>::operator = (matrix&& rhs){
        this->swap(rhs);
        return *this;
    }
#endif

    template<typename T>
    template<class M> 
    size_t matrix<T>::inc(const M& a){ 
        return 1; 
    }

    template<typename T>
    template<class M> 
    size_t matrix<T>::rows(const M& a){ 
        return a.num_rows(); 
    }

    template<typename T>
    template<class M> 
    size_t matrix<T>::cols(const M& a){ 
        return a.num_cols(); 
    }

    template<typename T>
    inline size_type matrix<T>::lda() const { 
        return this->core->spec.dim.y; 
    }

    template<typename T>
    inline size_type matrix<T>::num_rows() const { 
        return this->core->spec.dim.y; 
    }

    template<typename T>
    inline size_type matrix<T>::num_cols() const {
        return this->core->spec.dim.x; 
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
        return (this->core->spec.dim == 0);    
    }

    template<typename T>
    inline void matrix<T>::swap(matrix& r){ 
        std::swap(this->core, r.core);
    }

    template<typename T>
    inline void matrix<T>::resize(size_type m, size_type n){
        ambient::numeric::resize(*this, m, n);
    }

    template<typename T>
    inline void matrix<T>::remove_rows(size_type i, size_type k){
        ambient::numeric::remove_rows(*this, i, k);
    }

    template<typename T>
    inline void matrix<T>::remove_cols(size_type j, size_type k){
        ambient::numeric::remove_cols(*this, j, k); 
    }

    template<typename T>
    inline matrix<T>& matrix<T>::locate(size_type i, size_type j){
        return *this;
    }

    template<typename T>
    inline const matrix<T>& matrix<T>::locate(size_type i, size_type j) const {
        return *this;
    }

    template<typename T>
    inline size_t matrix<T>::addr(size_type i, size_type j) const {
        return (i + j*this->lda());
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
        ambient::model.touch(core);
        assert(ambient::model.common(r));
        ambient::sync(); if(!core->current->valid()) ambient::controller.calloc(*core->current); return ((value_type*)*core->current)[ j*core->spec.dim.y + i ];
    }

    template<typename T>
    inline const value_type& matrix<T>::operator() (size_type i, size_type j) const {
        ambient::model.touch(core);
        assert(ambient::model.common(r));
        ambient::sync(); if(!core->current->valid()) ambient::controller.calloc(*core->current); return ((value_type*)*core->current)[ j*core->spec.dim.y + i ];
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
