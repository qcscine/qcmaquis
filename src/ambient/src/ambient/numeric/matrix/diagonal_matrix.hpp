#ifndef AMBIENT_NUMERIC_DIAGONAL_MATRIX_HPP
#define AMBIENT_NUMERIC_DIAGONAL_MATRIX_HPP
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
    inline void diagonal_matrix<T>::remove_rows(size_t i, size_t m){
       this->data_.remove_rows(i, m);
    }

    template<typename T>  
    inline void diagonal_matrix<T>::remove_cols(size_t j, size_t n){
        this->data_.remove_rows(j, n);
    }

    template<typename T>
    inline diagonal_matrix<T>& diagonal_matrix<T>::locate(size_t i, size_t j){
        return *this;
    }

    template<typename T>
    inline const diagonal_matrix<T>& diagonal_matrix<T>::locate(size_t i, size_t j) const {
        return *this;
    }

    template<typename T>
    inline size_t diagonal_matrix<T>::addr(size_t i, size_t j) const {
        return i;
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

    // {{{ diagonal_matrix free functions
    template<typename T>
    inline size_type num_rows(const diagonal_matrix<T>& a){
        return a.num_rows();
    }

    template<typename T>
    inline size_type num_cols(const diagonal_matrix<T>& a){
        return a.num_cols();
    }

    template<typename T>
    inline diagonal_matrix< std::complex<T> > exp(const diagonal_matrix<T>& a, const std::complex<T>& alfa){
        assert(false); printf("ERROR: NOT TESTED (EXP)\n");
        diagonal_matrix< std::complex<T> > e(num_rows(a), num_rows(a));
        //kernels::exp_diagonal_rc<T>::spawn(e, a, alfa);
        return e;
    }

    template<typename T>
    inline diagonal_matrix<T> exp(diagonal_matrix<T> a, const T& alfa = 1.){
        exp_inplace(a, alfa);
        return a;
    }

    template<typename T>
    inline void exp_inplace(diagonal_matrix<T>& a, const T& alfa = 1.){
        kernels::exp_diagonal<T>::spawn<complexity::N>(a, alfa);
    }

    template<typename T>
    inline diagonal_matrix<T> sqrt(diagonal_matrix<T> a){
        sqrt_inplace(a);
        return a;
    }

    template<typename T>
    inline void sqrt_inplace(diagonal_matrix<T>& a){
        kernels::sqrt_diagonal<T>::spawn<complexity::N>(a);
    }

    template<typename T>
    inline std::ostream& operator<<(std::ostream& os, const diagonal_matrix<T>& a){
       os << a.data_ ; 
       return os;
    }

    template<typename T>
    inline void remove_rows(diagonal_matrix<T>& a, size_t i, size_t m = 1){
        a.remove_rows(i, m);
    }

    template<typename T>
    inline void remove_cols(diagonal_matrix<T>& a, size_t j, size_t n = 1){
        a.remove_cols(j, n);
    }

    template<typename T>
    inline void resize(diagonal_matrix<T>& a, size_t m, size_t n){
        a.resize(m, n);
    }
    // }}}

    #undef value_type
    #undef size_type

} } // namespace ambient::numeric

#endif
