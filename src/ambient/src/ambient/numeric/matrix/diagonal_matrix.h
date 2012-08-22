#ifndef __AMBIENT_NUMERIC_DIAGONAL_MATRIX_H
#define __AMBIENT_NUMERIC_DIAGONAL_MATRIX_H

namespace ambient { namespace numeric {
   
    template<typename T>
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
        
        inline diagonal_matrix(size_t rows, const value_type& init = value_type());
        inline size_type num_rows() const;
        inline size_type num_cols() const;
        inline const value_type& operator[](size_t i) const;
        inline value_type& operator[](size_t i); 
        inline const value_type& operator()(size_t i, size_t j) const;
        inline value_type& operator()(size_t i, size_t j);
        inline void remove_rows(size_t i, size_t k = 1);
        inline void remove_cols(size_t j, size_t k = 1);
        inline void resize(size_t rows, size_t cols);
        template< class T1 > 
        friend std::ostream & operator<<(std::ostream& os, const diagonal_matrix<T1>& m);
        inline const container& get_data() const; 
        inline container& get_data();    
        inline size_type size() const;
        inline void exp(const T& alfa = 1.);
   private:
        container data_;
    };

} } // namespace ambient::numeric

#endif
