#include "types/p_dense_matrix/p_diagonal_matrix.h" 

namespace maquis { namespace types {

    #define size_type   typename p_diagonal_matrix<T>::size_type
    #define value_type  typename p_diagonal_matrix<T>::value_type

    template<typename T>
    inline p_diagonal_matrix<T>::p_diagonal_matrix(size_t rows, const value_type& init)
    : data_(rows,1,init)
    {
    }

    template<typename T>
    inline size_type p_diagonal_matrix<T>::num_rows() const {
        return this->data_.num_rows();
    }

    template<typename T>
    inline size_type p_diagonal_matrix<T>::num_cols() const {
        return this->num_rows();
    }

    template<typename T>
    inline const value_type& p_diagonal_matrix<T>::operator[](size_t i) const {
        return this->data_(i,0);
    }

    template<typename T>
    inline value_type& p_diagonal_matrix<T>::operator[](size_t i){
        return this->data_(i,0);
    }

    template<typename T>
    inline const value_type& p_diagonal_matrix<T>::operator()(size_t i, size_t j) const {
        assert(i == j);
        return this->data_(i,0);
    }

    template<typename T>
    inline value_type& p_diagonal_matrix<T>:: operator()(size_t i, size_t j){
        assert(i == j);
        return this->data_(i,0);
    }

    template<typename T>
    inline size_type p_diagonal_matrix<T>::size() const {
        return this->data_.num_rows();
    }

    template<typename T>
    inline void p_diagonal_matrix<T>::remove_rows(size_t i, size_t k){
       this->data_.remove_rows(i, k);
    }

    template<typename T>  
    inline void p_diagonal_matrix<T>::remove_cols(size_t j, size_t k){
        this->data_.remove_rows(j, k);
    }

    template<typename T> 
    inline void p_diagonal_matrix<T>::resize(size_t rows, size_t cols){
        assert(rows == cols);
        this->data_.resize(rows, 1);
    }

    template<typename T>
    inline const typename p_diagonal_matrix<T>::container& p_diagonal_matrix<T>::get_data() const {
        return this->data_;                          
    }

    template<typename T>
    inline typename p_diagonal_matrix<T>::container& p_diagonal_matrix<T>::get_data(){
        return this->data_;                          
    }

    template<typename T>
    inline void p_diagonal_matrix<T>::exp(const T& alfa){
        assert(false); printf("NOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO <- EXP DIAG");
        ambient::push< ambient::exp_diagonal<T> >(*this, alfa);
    }

    template<typename T>
    inline void p_diagonal_matrix<T>::sqrt(){ 
        assert(false); printf("NOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO <- SQRT DIAG");
        ambient::push< ambient::sqrt_diagonal<T> >(*this);
    }

    // {{{ p_diagonal_matrix free functions
    template<typename T>
    inline size_type num_rows(const p_diagonal_matrix<T>& m){
        return m.num_rows();
    }

    template<typename T>
    inline size_type num_cols(const p_diagonal_matrix<T>& m){
        return m.num_cols();
    }

    template<typename T>
    inline p_diagonal_matrix<T> exp(const p_diagonal_matrix<T>& m, const T& alfa = 1.){
        p_diagonal_matrix<T> result(m);
        result.exp(alfa);
        return result;
    }

    template<typename T>
    inline p_diagonal_matrix< std::complex<T> > exp(const p_diagonal_matrix<T>& m, const std::complex<T>& alfa){
        assert(false); printf("NOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO <- EXP");
        assert(false);
        p_diagonal_matrix< std::complex<T> > e(num_rows(m), num_rows(m));
        ambient::push< ambient::exp_diagonal_rc<T> >(e, m, alfa);
        return e;
    }

    template<typename T>
    inline p_diagonal_matrix<T> sqrt(const p_diagonal_matrix<T>& m){
        p_diagonal_matrix<T> result(m);
        result.sqrt();
        return result;
    }

    template<typename T>
    inline std::ostream& operator<<(std::ostream& os, const p_diagonal_matrix<T>& m){
       os << m.data_ ; 
       return os;
    }

    template<typename T>
    inline void remove_rows(p_diagonal_matrix<T>& m, size_t k, size_t n = 1){
        m.remove_rows(k, n);
    }

    template<typename T>
    inline void remove_cols(p_diagonal_matrix<T>& m, size_t k, size_t n = 1){
        m.remove_cols(k, n);
    }

    template<typename T>
    inline void resize(p_diagonal_matrix<T>& m, size_t rows, size_t cols){
        m.resize(rows, cols);
    }
    // }}}

    #undef value_type
    #undef size_type

} } // namespace maquis::types

