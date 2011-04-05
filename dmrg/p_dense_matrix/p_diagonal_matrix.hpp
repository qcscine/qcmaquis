#include "p_dense_matrix/p_diagonal_matrix.h" 


namespace blas {

    template<typename T>
    p_diagonal_matrix<T>::p_diagonal_matrix(size_t rows, T const & init):data_(rows,1)  
    {
    }

    template<typename T>
    std::size_t p_diagonal_matrix<T>::num_rows() const 
    {
       return this->data_.num_rows();
    }

    template<typename T>
    std::size_t p_diagonal_matrix<T>::num_cols() const
    {
        return this->num_rows();
    }
    
    template<typename T>
    T const & p_diagonal_matrix<T>::operator[](std::size_t i) const 
    {
        return this->data_(i,0);
    }

    template<typename T>
    T & p_diagonal_matrix<T>::operator[](std::size_t i)
    {
        return this->data_(i,0);
    }
    
    template<typename T>
    T const & p_diagonal_matrix<T>::operator()(std::size_t i, std::size_t j) const
    {
        assert(i == j);
        return this->data_(i,0);
    }

    template<typename T>
    T & p_diagonal_matrix<T>:: operator()(std::size_t i, std::size_t j)
    {
        assert(i == j);
        return this->data_(i,0);
    }

    template<typename T>
    std::pair<typename p_diagonal_matrix<T>::element_iterator, typename p_diagonal_matrix<T>::element_iterator> p_diagonal_matrix<T>::elements()
    {
        int n = this->data_.num_rows(); 
        return std::make_pair(this->data_(0,0), this->data_(n,0));
    }

    template<typename T>
    std::pair<typename p_diagonal_matrix<T>::const_element_iterator, typename p_diagonal_matrix<T>::const_element_iterator> p_diagonal_matrix<T>::elements() const
    {
        int n = this->data_.num_rows(); 
        return std::make_pair(this->data_(0,0), this->data_(n,0));
    }
        


    template<typename T>
    void p_diagonal_matrix<T>::remove_rows(std::size_t k, std::size_t n)
    {
       this->data_.remove_rows(k, n);
    }

    template<typename T>  
    void p_diagonal_matrix<T>::remove_cols(std::size_t k, std::size_t n)
    {
        this->data_.remove_rows(k, n);
    }

    template<typename T> 
    void p_diagonal_matrix<T>::resize(std::size_t rows, std::size_t cols, T v)
    {
        assert(rows == cols);
        this->data_.resize(rows, 1);
    }

    template<typename T>
    ambient::p_dense_matrix<T> const  & p_diagonal_matrix<T>::get_data() const 
    {
        return this->data_;                          
    }                                             

    template<typename T>
    ambient::p_dense_matrix<T> & p_diagonal_matrix<T>::get_data() 
    {
        return this->data_;                          
    }

    template<typename T, class Matrix>
    void gemm(Matrix const & m1, p_diagonal_matrix<T> const & m2, Matrix & m3)
    {
        assert(num_cols(m1) == num_rows(m2));
        m3.resize(num_rows(m1), num_cols(m2));
	ambient::push(ambient::gemm_rhs_diagonal_l_kernel, ambient::gemm_rhs_diagonal_c_kernel, m1, m2.get_data() ,m3);
        /*
        for (std::size_t i = 0; i < num_rows(m1); ++i)
            for (std::size_t j = 0; j < num_columns(m2); ++j)
                m3(i,j) = m1(i,j) * m2(j,j);
        */
    }
    
//to do wrapper scala first
    template<typename T, class Matrix>
    void gemm( const p_diagonal_matrix<T> & m1, Matrix const & m2, Matrix & m3)
    {
        assert(num_columns(m1) == num_rows(m2));
        m3.resize(num_rows(m1), num_columns(m2));
	//ambient::push(ambient::gemm_rhs_diagonal_l_kernel,ambient::gemm_rhs_diagonal_c_kernel, m2, m1.get_data() ,m3)
//	ambient::push(ambient::gemm_rhs_diagonal_l_kernel,ambient::gemm_rhs_diagonal_c_kernel, m2, m3)
/*
        for (std::size_t i = 0; i < num_rows(m1); ++i)
            for (std::size_t j = 0; j < num_columns(m2); ++j)
                m3(i,j) = m1(i,i) * m2(i,j);
  */  
    }

    template<typename T>
    typename p_diagonal_matrix<T>::size_type num_rows(p_diagonal_matrix<T> const & m)
    {
        return m.num_rows();
    }
    
    template<typename T>
    typename p_diagonal_matrix<T>::size_type num_cols(p_diagonal_matrix<T> const & m)
    {
        return m.num_cols();
    }
 

    template<typename T>
    void sqrt(p_diagonal_matrix<T> & m)
    { 
        ambient::push(ambient::sqrt_diagonal_l_kernel, ambient::sqrt_diagonal_c_kernel, m.get_data());
    }
   
    template<typename T>
    std::ostream& operator<<(std::ostream& os, p_diagonal_matrix<T> const & m)
    {
       os << m.data_ ; 
       return os;
    }
 
    template<typename T>
    void remove_rows(p_diagonal_matrix<T> & m, std::size_t k, std::size_t n = 1)
    {
        m.remove_rows(k, n);
    }
    
    template<typename T>
    void remove_cols(p_diagonal_matrix<T> & m, std::size_t k, std::size_t n = 1)
    {
        m.remove_cols(k, n);
    }
   

    template<typename T>
    void resize(p_diagonal_matrix<T> & m, std::size_t rows, std::size_t cols, T v = T())
    {
        m.resize(rows, cols, v);
    }
    
    template<typename T>
    std::pair<typename p_diagonal_matrix<T>::element_iterator, typename p_diagonal_matrix<T>::element_iterator>
    elements(p_diagonal_matrix<T> & m)
    {
        return m.elements();
    }

    template<typename T>
    std::pair<typename p_diagonal_matrix<T>::const_element_iterator, typename p_diagonal_matrix<T>::const_element_iterator>
    elements(p_diagonal_matrix<T> const & m)
    {
        return m.elements();
    }
}

