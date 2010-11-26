#ifndef DIAGONAL_MATRIX_H
#define DIAGONAL_MATRIX_H

namespace blas {

    template<typename T>
    class diagonal_matrix
    {
    public:
        typedef T                       value_type;
        typedef T&                      reference;
        typedef T const&                const_reference;
        typedef std::size_t             size_type;
        typedef std::ptrdiff_t          difference_type;
        
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
        
    private:
        std::vector<T> data_;
    };
    
    template<typename T, class Matrix>
    void gemm(Matrix const & m1, diagonal_matrix<T> const & m2, Matrix & m3)
    {
        assert(num_columns(m1) == num_rows(m2));
        resize(m3, num_rows(m1), num_columns(m2));
        for (std::size_t i = 0; i < num_rows(m1); ++i)
            for (std::size_t j = 0; j < num_rows(m2); ++j)
                m3(i,j) = m1(i,j) * m2(j,j);
    }
    
    template<typename T, class Matrix>
    void gemm(diagonal_matrix<T> const & m1, Matrix const & m2, Matrix & m3)
    {
        assert(num_columns(m1) == num_rows(m2));
        resize(m3, num_rows(m1), num_columns(m2));
        for (std::size_t i = 0; i < num_rows(m1); ++i)
            for (std::size_t j = 0; j < num_rows(m2); ++j)
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
}

#endif
