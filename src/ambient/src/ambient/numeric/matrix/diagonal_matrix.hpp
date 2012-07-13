#ifndef __AMBIENT_NUMERIC_DIAGONAL_MATRIX_HPP
#define __AMBIENT_NUMERIC_DIAGONAL_MATRIX_HPP
#include "ambient/numeric/matrix/diagonal_matrix.h" 

namespace ambient { namespace numeric {

    #define size_type   typename diagonal_matrix<T>::size_type
    #define value_type  typename diagonal_matrix<T>::value_type

    template<typename T>
    inline diagonal_matrix<T>::diagonal_matrix(size_t rows, const value_type& init)
    : data_(rows,1,init)
    {
    }

    template<typename T>
    inline size_type diagonal_matrix<T>::num_rows() const {
        return this->data_.num_rows();
    }

    template<typename T>
    inline size_type diagonal_matrix<T>::num_cols() const {
        return this->num_rows();
    }

    template<typename T>
    inline const value_type& diagonal_matrix<T>::operator[](size_t i) const {
        return this->data_(i,0);
    }

    template<typename T>
    inline value_type& diagonal_matrix<T>::operator[](size_t i){
        return this->data_(i,0);
    }

    template<typename T>
    inline const value_type& diagonal_matrix<T>::operator()(size_t i, size_t j) const {
        assert(i == j);
        return this->data_(i,0);
    }

    template<typename T>
    inline value_type& diagonal_matrix<T>:: operator()(size_t i, size_t j){
        assert(i == j);
        return this->data_(i,0);
    }

    template<typename T>
    inline size_type diagonal_matrix<T>::size() const {
        return this->data_.num_rows();
    }

    template<typename T>
    inline void diagonal_matrix<T>::remove_rows(size_t i, size_t k){
       this->data_.remove_rows(i, k);
    }

    template<typename T>  
    inline void diagonal_matrix<T>::remove_cols(size_t j, size_t k){
        this->data_.remove_rows(j, k);
    }

    template<typename T> 
    inline void diagonal_matrix<T>::resize(size_t rows, size_t cols){
        assert(rows == cols);
        this->data_.resize(rows, 1);
    }

    template<typename T>
    inline const typename diagonal_matrix<T>::container& diagonal_matrix<T>::get_data() const {
        return this->data_;                          
    }

    template<typename T>
    inline typename diagonal_matrix<T>::container& diagonal_matrix<T>::get_data(){
        return this->data_;                          
    }

    template<typename T>
    inline void diagonal_matrix<T>::exp(const T& alfa){
        assert(false); printf("ERROR: NOT TESTED (EXP DIAG)\n");
        //ambient::push< kernels::exp_diagonal<T> >(*this, alfa);
    }

    // {{{ diagonal_matrix free functions
    template<typename T>
    inline size_type num_rows(const diagonal_matrix<T>& m){
        return m.num_rows();
    }

    template<typename T>
    inline size_type num_cols(const diagonal_matrix<T>& m){
        return m.num_cols();
    }

    template<typename T>
    inline diagonal_matrix<T> exp(diagonal_matrix<T> m, const T& alfa = 1.){
        m.exp(alfa);
        return m;
    }

    template<typename T>
    inline diagonal_matrix< std::complex<T> > exp(const diagonal_matrix<T>& m, const std::complex<T>& alfa){
        assert(false); printf("ERROR: NOT TESTED (EXP)\n");
        diagonal_matrix< std::complex<T> > e(num_rows(m), num_rows(m));
        //ambient::push< kernels::exp_diagonal_rc<T> >(e, m, alfa);
        return e;
    }

    template<typename T>
    inline diagonal_matrix<T> sqrt(diagonal_matrix<T> m){
        sqrt_inplace(m);
        return m;
    }

    template<typename T>
    inline void sqrt_inplace(diagonal_matrix<T>& m){
        assert(false); printf("ERROR: NOT TESTED (SQRT DIAG)\n");
        //ambient::push< kernels::sqrt_diagonal<T> >(m);
    }

    template<typename T>
    inline std::ostream& operator<<(std::ostream& os, const diagonal_matrix<T>& m){
       os << m.data_ ; 
       return os;
    }

    template<typename T>
    inline void remove_rows(diagonal_matrix<T>& m, size_t k, size_t n = 1){
        m.remove_rows(k, n);
    }

    template<typename T>
    inline void remove_cols(diagonal_matrix<T>& m, size_t k, size_t n = 1){
        m.remove_cols(k, n);
    }

    template<typename T>
    inline void resize(diagonal_matrix<T>& m, size_t rows, size_t cols){
        m.resize(rows, cols);
    }
    // }}}

    #undef value_type
    #undef size_type

} } // namespace ambient::numeric

#endif
