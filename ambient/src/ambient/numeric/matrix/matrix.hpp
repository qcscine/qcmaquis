#ifndef AMBIENT_NUMERIC_MATRIX_HPP
#define AMBIENT_NUMERIC_MATRIX_HPP

#include "ambient/numeric/matrix/matrix.h"
#include "ambient/numeric/matrix/matrix_algorithms.hpp"
#ifdef AMBIENT_EXPERIMENTAL
#include "ambient/numeric/matrix/matrix_algorithms_experimental.hpp"
#endif


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
    #define size_type   typename matrix<T,A>::size_type
    #define value_type  typename matrix<T,A>::value_type
    #define scalar_type typename matrix<T,A>::scalar_type

    template<typename T, class A>
    inline void* matrix<T,A>::operator new (size_t size){
        return ambient::pool.malloc<matrix<T,A> >();
    }

    template<typename T, class A>
    inline void* matrix<T,A>::operator new (size_t size, void* placement){
        return placement; 
    }

    template<typename T, class A>
    inline void matrix<T,A>::operator delete (void* ptr){
        ambient::pool.free<matrix<T,A> >(ptr);
    }

    template <typename T, class A>
    inline matrix<T,A>::~matrix(){
        if(this->core->weak()) delete this->core;
        else ambient::destroy(this->core);
    }

    template <typename T, class A>
    inline matrix<T,A>::matrix(const ptr& p, size_t r) 
    : core(p), ref(r)
    {
    }

    template <typename T, class A>
    inline matrix<T,A>::matrix(){ 
        this->core = new I(ambient::dim2(0,0), sizeof(T)); 
    }

    template <typename T, class A>
    inline matrix<T,A>::matrix(size_type rows, size_type cols, value_type init_value){
        this->core = new I(ambient::dim2(cols, rows), sizeof(T));
        fill_value(*this, init_value);
    }

    template <typename T, class A>
    inline matrix<T,A>::matrix(const matrix& a){
        this->core = new I(a.core->spec.dim, sizeof(T));
        ambient::fuse(a.core, this->core);
    }
    
    template <typename T, class A>
    matrix<T,A>& matrix<T,A>::operator = (const matrix& rhs){
        //if(rhs.core->weak() && this->core->weak()){
        //    this->core->spec.dim = rhs.core->spec.dim;
        //    return *this;
        //}   // commented due to insignificance
        matrix c(rhs);
        this->swap(c);
        return *this;
    }

    template<typename T, class A>
    template<class M> 
    size_t matrix<T,A>::inc(const M& a){ 
        return 1; 
    }

    template<typename T, class A>
    template<class M> 
    size_t matrix<T,A>::rows(const M& a){ 
        return a.num_rows(); 
    }

    template<typename T, class A>
    template<class M> 
    size_t matrix<T,A>::cols(const M& a){ 
        return a.num_cols(); 
    }

    template<typename T, class A>
    inline size_type matrix<T,A>::lda() const { 
        return this->core->spec.dim.y; 
    }

    template<typename T, class A>
    inline size_type matrix<T,A>::num_rows() const { 
        return this->core->spec.dim.y; 
    }

    template<typename T, class A>
    inline size_type matrix<T,A>::num_cols() const {
        return this->core->spec.dim.x; 
    }

    template<typename T, class A>
    inline scalar_type matrix<T,A>::trace() const { 
        return trace(*this);           
    }

    template<typename T, class A>
    inline void matrix<T,A>::transpose(){ 
        transpose_inplace(*this);      
    }

    template<typename T, class A>
    inline void matrix<T,A>::conj(){ 
        conj_inplace(*this);           
    }

    template<typename T, class A>
    inline bool matrix<T,A>::empty() const { 
        return (this->core->spec.dim == 0);    
    }

    template<typename T, class A>
    inline void matrix<T,A>::swap(matrix& r){ 
        std::swap(this->core, r.core);
    }

    template<typename T, class A>
    inline void matrix<T,A>::resize(size_type m, size_type n){
        ambient::numeric::resize(*this, m, n);
    }

    template<typename T, class A>
    inline void matrix<T,A>::remove_rows(size_type i, size_type k){
        ambient::numeric::remove_rows(*this, i, k);
    }

    template<typename T, class A>
    inline void matrix<T,A>::remove_cols(size_type j, size_type k){
        ambient::numeric::remove_cols(*this, j, k); 
    }

    template<typename T, class A>
    inline matrix<T,A>& matrix<T,A>::locate(size_type i, size_type j){
        return *this;
    }

    template<typename T, class A>
    inline const matrix<T,A>& matrix<T,A>::locate(size_type i, size_type j) const {
        return *this;
    }

    template<typename T, class A>
    inline size_t matrix<T,A>::addr(size_type i, size_type j) const {
        return (i + j*this->lda());
    }

    template<typename T, class A>
    inline matrix<T,A>& matrix<T,A>::operator += (const matrix& rhs){
        add_inplace(*this, rhs);
        return *this;
    }

    template<typename T, class A>
    inline matrix<T,A>& matrix<T,A>::operator -= (const matrix& rhs){
        sub_inplace(*this, rhs);
        return *this;
    }

    template<typename T, class A>
    template <typename T2> 
    inline matrix<T,A>& matrix<T,A>::operator *= (const T2& t){
        mul_inplace(*this, t);
        return *this;
    }

    template<typename T, class A>
    template <typename T2> 
    inline matrix<T,A>& matrix<T,A>::operator /= (const T2& t){
        div_inplace(*this, t);
        return *this;
    }

    template<typename T, class A>
    inline value_type& matrix<T,A>::operator() (size_type i, size_type j){
        return ((T*)ambient::serial(*this))[ j*lda() + i ];
    }

    template<typename T, class A>
    inline const value_type& matrix<T,A>::operator() (size_type i, size_type j) const {
        return ((T*)ambient::serial(*this))[ j*lda() + i ];
    }

    template<typename T, class A>
    const char* matrix<T,A>::code(){ 
        return "N"; 
    }  
    #undef size_type
    #undef value_type
    #undef scalar_type
    // }}}

} }

#endif
