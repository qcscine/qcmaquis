#include "types/p_dense_matrix/p_dense_matrix.h"
#include "types/p_dense_matrix/algorithms.hpp"

namespace maquis { namespace types {

    #define size_type typename p_dense_matrix_impl<T>::size_type

    template <typename T>
    p_dense_matrix_impl<T>::~p_dense_matrix_impl(){ } // #destructor // 

    template <typename T>
    p_dense_matrix_impl<T>::p_dense_matrix_impl(size_type rows, size_type cols = 0, T init_value = T() ){
        this->cols = cols;
        this->rows = rows;
        this->init_value = init_value;
        this->pt_set_dim(cols, rows);
        this->pt_set_init(ambient::value_i<T>);
    }

    template <typename T>
    p_dense_matrix_impl<T>::p_dense_matrix_impl(const p_dense_matrix_impl& m)
    : ambient::parallel_t< p_dense_matrix_impl<T> >(m)
    {
        this->cols = m.num_cols();
        this->rows = m.num_rows();
    }

    template <typename T>
    inline bool p_dense_matrix_impl<T>::empty() const {
        return (this->rows == 0 || this->cols == 0); 
    }

    template <typename T>
    inline size_type p_dense_matrix_impl<T>::num_rows() const {
        return this->rows;   
    }

    template <typename T>
    inline size_type p_dense_matrix_impl<T>::num_cols() const {
        return this->cols;   
    }

    template <typename T>
    void p_dense_matrix_impl<T>::clear(){
        algorithms::clear(*this);
        this->rows = this->cols = 0;
    }

    template <typename T>
    void p_dense_matrix_impl<T>::resize(size_type rows, size_type cols){
        assert(rows > 0); assert(cols > 0);
        algorithms::resize(*this, rows, cols);
        this->rows = rows; 
        this->cols = cols;
    }

    template <typename T>
    void p_dense_matrix_impl<T>::remove_rows(size_type i, size_type k = 1){
        assert( i+k <= this->rows );
        algorithms::remove_rows(*this, i, k);
        this->rows -= k; // no storage shrink
    }

    template <typename T>
    void p_dense_matrix_impl<T>::remove_cols(size_type j, size_type k = 1){
        assert( j+k <= this->cols );
        algorithms::remove_cols(*this, j, k);
        this->cols -= k; // no storage shrink
    }

    template <typename T>
    void p_dense_matrix_impl<T>::conjugate(){ 
        algorithms::inplace_conjugate(*this);
    }

    template <typename T>
    void p_dense_matrix_impl<T>::transpose(){ 
        algorithms::inplace_transpose(*this);
    }

    template <typename T>
    T& p_dense_matrix_impl<T>::get(size_type i, size_type j){
        assert(i < this->rows); assert(j < this->cols);
        return this->pt_fetch(i / this->pt_mem_dim().y,  // blocked_i
                              j / this->pt_mem_dim().x,  // blocked_j 
                              i % this->pt_mem_dim().y,  // element_i 
                              j % this->pt_mem_dim().x); // element_j
    }

    template <typename T>
    T p_dense_matrix_impl<T>::trace() {
        return algorithms::trace(*this);
    }

    template <typename T>
    void p_dense_matrix_impl<T>::add(const p_dense_matrix_impl& rhs){
        algorithms::add_inplace(*this, rhs);
    }

    template <typename T>
    void p_dense_matrix_impl<T>::sub(const p_dense_matrix_impl& rhs){
        algorithms::sub_inplace(*this, rhs);
    }

    template <typename T>
    void p_dense_matrix_impl<T>::mul(const p_dense_matrix_impl<T>& rhs){
        algorithms::gemm_inplace(*this, rhs);
    }

    template <typename T>
    template <typename T2>
    void p_dense_matrix_impl<T>::mul(const T2& t){
        algorithms::scale_inplace(*this, t);
    }

    template <typename T>
    void p_dense_matrix_impl<T>::cpy(const p_dense_matrix_impl<T>& rhs){
        this->resize(rhs.num_rows(), rhs.num_cols());
        algorithms::cpy(*this, rhs);
    }

    #undef size_type

} } // namespace maquis::types

