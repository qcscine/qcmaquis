#ifndef __AMBIENT_NUMERIC_MATRIX_HPP__
#define __AMBIENT_NUMERIC_MATRIX_HPP__
#include "ambient/numeric/matrix/matrix.h"
#include "ambient/numeric/matrix/algorithms.hpp"

namespace ambient { namespace numeric {

    #define size_type   typename matrix_impl<T>::size_type
    #define value_type  typename matrix_impl<T>::value_type
    #define scalar_type typename matrix_impl<T>::scalar_type

    template <typename T>
    inline matrix_impl<T>::matrix_impl()
    : ambient::iteratable<history>(ambient::dim2(0,0)), references(0)
    { // be cautious (implicit)
    }

    template <typename T>
    inline matrix_impl<T>::matrix_impl(size_type rows, size_type cols)
    : ambient::iteratable<history>(ambient::dim2(cols, rows)), references(0)
    {
    }

    template <typename T>
    inline matrix_impl<T>::matrix_impl(const matrix_impl& m)
    : ambient::iteratable<history>(m.spec.dim), references(0)
    {
    }

    template <typename T>
    inline bool matrix_impl<T>::empty() const {
        return (this->spec.dim.x == 0 || this->spec.dim.y == 0); 
    }

    template <typename T>
    inline size_type matrix_impl<T>::num_rows() const {
        return this->spec.dim.y;   
    }

    template <typename T>
    inline size_type matrix_impl<T>::num_cols() const {
        return this->spec.dim.x;   
    }

    template <typename T>
    inline bool matrix_impl<T>::atomic() const {
        return ambient::model.is_atomic(this);
    }

    template <typename T>
    inline void matrix_impl<T>::resize(matrix_impl& r, size_type rows, size_type cols){
        algorithms::resize(r, rows, cols, *this, this->spec.dim.y, this->spec.dim.x);
    }

    template <typename T>
    inline void matrix_impl<T>::remove_rows(size_type i, size_type k){
        assert( i+k <= this->spec.dim.y );
        algorithms::remove_rows(*this, i, k);
    }

    template <typename T>
    inline void matrix_impl<T>::remove_cols(size_type j, size_type k){
        assert( j+k <= this->spec.dim.x );
        algorithms::remove_cols(*this, j, k);
    }

    template <typename T>
    inline void matrix_impl<T>::fill_identity(){ 
        algorithms::fill_identity(*this);
    }

    template <typename T>
    inline void matrix_impl<T>::fill_random(){ 
        algorithms::fill_random(*this);
    }

    template <typename T>
    inline void matrix_impl<T>::fill_value(value_type v){ 
        algorithms::fill_value(*this, v);
    }

    template <typename T>
    inline void matrix_impl<T>::conj(){ 
        algorithms::conj_inplace(*this);
    }

    template <typename T>
    inline void matrix_impl<T>::transpose(){ 
        algorithms::transpose_inplace(*this);
    }

    template <typename T>
    inline value_type& matrix_impl<T>::get(size_type i, size_type j){
        printf("Get element explicit PLAYOUT!\n");
        ambient::playout();
        return ((value_type*)ambient::controller.ufetch_block(*this->current, j/this->spec.block.x, i/this->spec.block.y))
               [ (j%this->spec.block.x)*this->spec.block.y + i%this->spec.block.y ];
    }

    template <typename T>
    inline scalar_type matrix_impl<T>::trace() const {
        return algorithms::trace(*this);
    }

    template <typename T>
    inline void matrix_impl<T>::add(const matrix_impl& rhs){
        algorithms::add_inplace(*this, rhs);
    }

    template <typename T>
    inline void matrix_impl<T>::sub(const matrix_impl& rhs){
        algorithms::sub_inplace(*this, rhs);
    }

    template <typename T>
    inline void matrix_impl<T>::mul(const diagonal_matrix<T>& rhs){
        algorithms::gemm_diag_inplace(*this, *rhs.get_data().impl);
    }

    template <typename T>
    inline void matrix_impl<T>::mul(const matrix_impl<T>& rhs){
        algorithms::gemm_inplace(*this, rhs);
    }

    template <typename T>
    template <typename T2>
    inline void matrix_impl<T>::mul(const T2& t){
        algorithms::scale_inplace(*this, t);
    }

    template <typename T>
    inline void matrix_impl<T>::copy(const matrix_impl<T>& m){
        algorithms::copy(*this, m);
    }

    #undef size_type
    #undef value_type
    #undef scalar_type

} } // namespace ambient::numeric

#endif
