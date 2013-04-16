#ifndef AMBIENT_NUMERIC_DIAGONAL_MATRIX_H
#define AMBIENT_NUMERIC_DIAGONAL_MATRIX_H

namespace ambient { namespace numeric {
   
    template<typename T, class Allocator>
    class matrix;

    template<typename T>
    class diagonal_matrix
    {
    public:
        typedef matrix<T> container;
        typedef typename container::difference_type difference_type;
        typedef typename container::value_type value_type;
        typedef typename container::scalar_type scalar_type;
        typedef typename container::real_type real_type;
        typedef typename container::size_type size_type;
        
        diagonal_matrix(size_t rows, const value_type& init = value_type());
        size_type num_rows() const;
        size_type num_cols() const;
        const value_type& operator[](size_t i) const;
        value_type& operator[](size_t i); 
        const value_type& operator()(size_t i, size_t j) const;
        value_type& operator()(size_t i, size_t j);
        void remove_rows(size_t i, size_t m = 1);
        void remove_cols(size_t j, size_t n = 1);
        diagonal_matrix& locate(size_t i, size_t j);
        const diagonal_matrix& locate(size_t i, size_t j) const;
        size_t addr(size_t i, size_t j) const;
        void resize(size_t m, size_t n);
        template< class T1 > 
        friend std::ostream & operator<<(std::ostream& os, const diagonal_matrix<T1>& a);
        const container& get_data() const; 
        container& get_data();    
        size_type size() const;
   private:
        container data_;
    };

} }

#endif
