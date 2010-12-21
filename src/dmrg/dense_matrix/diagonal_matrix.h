#ifndef DIAGONAL_MATRIX_H
#define DIAGONAL_MATRIX_H

namespace blas {

    template<class FullMatrixClass>
    struct associated_diagonal_matrix { };
    
    template<typename T>
    class diagonal_matrix
    {
    public:
        typedef T                       value_type;
        typedef T&                      reference;
        typedef T const&                const_reference;
        typedef std::size_t             size_type;
        typedef std::ptrdiff_t          difference_type;
        
        typedef typename std::vector<T>::iterator element_iterator;
        typedef typename std::vector<T>::const_iterator const_element_iterator;
        
        template<class Vector>
        diagonal_matrix(Vector const & init)
        : data_(init.begin(), init.end()) { }

        diagonal_matrix(std::size_t size, T const & init = T())
        : data_(size, init) { }
        
        std::size_t num_rows() const { return data_.size(); }
        std::size_t num_columns() const { return data_.size(); }
        
        T const & operator[](std::size_t i) const { return data_[i]; }
        T & operator[](std::size_t i) { return data_[i]; }
        
        T const & operator()(std::size_t i, std::size_t j) const
        {
            assert(i == j);
            return data_[i];
        }
        T & operator()(std::size_t i, std::size_t j)
        {
            assert(i == j);
            return data_[i];
        }
        
        std::pair<element_iterator, element_iterator> elements()
        {
            return std::make_pair(data_.begin(), data_.end());
        }
        
        std::pair<const_element_iterator, const_element_iterator> elements() const
        {
            return std::make_pair(data_.begin(), data_.end());
        }
        
        void remove_rows(std::size_t k, std::size_t n = 1)
        {
            data_.erase(data_.begin(), data_.begin()+n);
        }
        
        void remove_columns(std::size_t k, std::size_t n)
        {
            remove_rows(k, n);
        }
        
        void resize(std::size_t r, std::size_t c, T v = T())
        {
            assert(r == c);
            data_.resize(r, v);
        }
        
    private:
        std::vector<T> data_;
    };
    
    template<typename T, class Matrix>
    void gemm(Matrix const & m1, diagonal_matrix<T> const & m2, Matrix & m3)
    {
        assert(num_columns(m1) == num_rows(m2));
        resize(m3, num_rows(m1), num_columns(m2));
        for (std::size_t i = 0; i < num_rows(m1); ++i)
            for (std::size_t j = 0; j < num_columns(m2); ++j)
                m3(i,j) = m1(i,j) * m2(j,j);
    }
    
    template<typename T, class Matrix>
    void gemm(diagonal_matrix<T> const & m1, Matrix const & m2, Matrix & m3)
    {
        assert(num_columns(m1) == num_rows(m2));
        resize(m3, num_rows(m1), num_columns(m2));
        for (std::size_t i = 0; i < num_rows(m1); ++i)
            for (std::size_t j = 0; j < num_columns(m2); ++j)
                m3(i,j) = m1(i,i) * m2(i,j);
    }
    
    template<typename T>
    typename diagonal_matrix<T>::size_type num_rows(diagonal_matrix<T> const & m)
    {
        return m.num_rows();
    }
    
    template<typename T>
    typename diagonal_matrix<T>::size_type num_columns(diagonal_matrix<T> const & m)
    {
        return m.num_columns();
    }
    
    template<typename T>
    diagonal_matrix<T> sqrt(diagonal_matrix<T> m)
    {
        std::transform(m.elements().first, m.elements().second, m.elements().first, utils::functor_sqrt());
        return m;
    }
    
    template<typename T>
    std::ostream& operator<<(std::ostream& os, diagonal_matrix<T> const & m)
    {
        std::copy(m.elements().first, m.elements().second, std::ostream_iterator<T>(os, " "));
        return os;
    }
    
    template<typename T>
    void remove_rows(diagonal_matrix<T> & m, std::size_t k, std::size_t n = 1)
    {
        m.remove_rows(k, n);
    }
    
    template<typename T>
    void remove_columns(diagonal_matrix<T> & m, std::size_t k, std::size_t n = 1)
    {
        m.remove_columns(k, n);
    }
    
    template<typename T>
    void resize(diagonal_matrix<T> & m, std::size_t r, std::size_t c, T v = T())
    {
        m.resize(r, c, v);
    }
}

#endif
