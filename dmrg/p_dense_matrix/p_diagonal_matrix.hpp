#include "p_dense_matrix/p_diagonal_matrix.h" 


namespace blas {

    template<typename T>
    p_diagonal_matrix<T>::p_diagonal_matrix(size_t rows, const T& init):data_(rows,1){ }

    template<typename T>
    size_t p_diagonal_matrix<T>::num_rows() const {
       return this->data_.num_rows();
    }
    template<typename T>
    size_t p_diagonal_matrix<T>::num_cols() const {
        return this->num_rows();
    }
    template<typename T>
    p_diagonal_matrix<T>& p_diagonal_matrix<T>::operator=(const p_diagonal_matrix<T>& rhs){
        this->data_ = rhs.get_data();    
    }
    template<typename T>
    const T& p_diagonal_matrix<T>::operator[](size_t i) const {
        return this->data_(i,0);
    }
    template<typename T>
    T & p_diagonal_matrix<T>::operator[](size_t i){
        return this->data_(i,0);
    }
    template<typename T>
    const T& p_diagonal_matrix<T>::operator()(size_t i, size_t j) const {
        assert(i == j);
        return this->data_(i,0);
    }
    template<typename T>
    T & p_diagonal_matrix<T>:: operator()(size_t i, size_t j){
        assert(i == j);
        return this->data_(i,0);
    }
    template<typename T>
    std::size_t p_diagonal_matrix<T>::size() const {
        return this->data_.num_rows();
    }
    template<typename T>
    void p_diagonal_matrix<T>::remove_rows(size_t i, size_t k){
       this->data_.remove_rows(i, k);
    }
    template<typename T>  
    void p_diagonal_matrix<T>::remove_cols(size_t j, size_t k){
        this->data_.remove_rows(j, k);
    }
    template<typename T> 
    void p_diagonal_matrix<T>::resize(size_t rows, size_t cols, T v){
        assert(rows == cols);
        this->data_.resize(rows, 1);
    }
    template<typename T>
    const typename p_diagonal_matrix<T>::container& p_diagonal_matrix<T>::get_data() const {
        return this->data_;                          
    }                                             
    template<typename T>
    typename p_diagonal_matrix<T>::container& p_diagonal_matrix<T>::get_data(){
        return this->data_;                          
    }


// free-functions (possibly should return copies or proxy-objects) <- to check

    template<typename T>
    typename p_diagonal_matrix<T>::size_type num_rows(const p_diagonal_matrix<T>& m)
    {
        return m.num_rows();
    }
    template<typename T>
    typename p_diagonal_matrix<T>::size_type num_cols(const p_diagonal_matrix<T>& m)
    {
        return m.num_cols();
    }
    template<typename T>
    void sqrt(p_diagonal_matrix<T> & m)
    { 
        ambient::push(ambient::sqrt_diagonal_l, ambient::sqrt_diagonal_c, m.get_data());
    }
    template<typename T>
    std::ostream& operator<<(std::ostream& os, const p_diagonal_matrix<T>& m)
    {
       os << m.data_ ; 
       return os;
    }
    template<typename T>
    void remove_rows(p_diagonal_matrix<T> & m, size_t k, size_t n = 1)
    {
        m.remove_rows(k, n);
    }
    template<typename T>
    void remove_cols(p_diagonal_matrix<T> & m, size_t k, size_t n = 1)
    {
        m.remove_cols(k, n);
    }
    template<typename T>
    void resize(p_diagonal_matrix<T> & m, size_t rows, size_t cols, T v = T())
    {
        m.resize(rows, cols, v);
    }
}

