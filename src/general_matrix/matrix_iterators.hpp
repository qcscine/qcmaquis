#ifndef BLAS_MATRIX_ITERATORS
#define BLAS_MATRIX_ITERATORS

#include <boost/iterator/iterator_facade.hpp>
#include <boost/static_assert.hpp>

template <typename Matrix, typename T>
class matrix_column_element_iterator;

// replace by strided_iterator<T>
// T* p;
// std::size_t stride;


template <typename Matrix, typename T>
class matrix_row_element_iterator : public boost::iterator_facade<
                                matrix_row_element_iterator<Matrix,T>,
                                T,
                                boost::random_access_traversal_tag,
                                T&,
                                typename Matrix::difference_type
                                >
{
    // iterates over matrix elements within the same column


    public:
        typedef T value_type;

        matrix_row_element_iterator(Matrix* matrix,typename Matrix::size_type row, typename Matrix::size_type col)
            : m(matrix), row_pos(row), col_pos(col)
        {
            // The value_type of the iterator must be the value_type of the matrix or const Matrix::value_type
            BOOST_STATIC_ASSERT( (boost::is_same<typename Matrix::value_type, T>::value
                                 || boost::is_same<const typename Matrix::value_type,T>::value) );
        }

        template<typename Matrix2, typename U>
        matrix_row_element_iterator(matrix_row_element_iterator<Matrix2,U> const& r)
            : m(r.m), row_pos(r.row_pos), col_pos(r.col_pos)
            {}

        template<typename Matrix2, typename U>
        explicit matrix_row_element_iterator(matrix_column_element_iterator<Matrix2,U> const& col_iter)
            : m(col_iter.m), row_pos(col_iter.row), col_pos(col_iter.col)
            {}

    private:
        friend class boost::iterator_core_access;
        template <typename,typename> friend class matrix_row_element_iterator;

        value_type& dereference() const
        { return m->operator()(row_pos,col_pos); }

        // iterators are equal if they point to the same column of the same matrix
        // WARNING: since the row position is not compared
        // two iterators can be equal although they point to different elements
        template <typename U>
        bool equal(matrix_row_element_iterator<Matrix,U> const& y) const
        {
            if(m == y.m && col_pos == y.col_pos)
                return true;
            else
                return false;
        }
        void increment()
        {
            ++col_pos;
        }
        void decrement()
        {
            --col_pos;
        }
        void advance(typename Matrix::difference_type n)
        {
            col_pos += n;
        }

        template <typename U>
        typename Matrix::difference_type distance_to(matrix_row_element_iterator<Matrix,U> const& z) const
        {
            return z.col_pos - col_pos;
        }



        typename Matrix::size_type row_pos;
        typename Matrix::size_type col_pos;
        Matrix* m;
        
};


// replace by T*
template <typename Matrix, typename T>
class matrix_column_element_iterator : public boost::iterator_facade<
                                   matrix_column_element_iterator<Matrix,T>,
                                   T,
                                   boost::random_access_traversal_tag,
                                   T&,
                                   typename Matrix::difference_type
                                   >
{
    // iterates over matrix elements within the same row
    

    public:
        typedef T value_type;

        matrix_column_element_iterator(Matrix* matrix,typename Matrix::size_type row, typename Matrix::size_type col)
            : m(matrix), row_pos(row), col_pos(col)
            {
                // The value_type of the iterator must be the value_type of the matrix or const Matrix::value_type
                BOOST_STATIC_ASSERT( (boost::is_same<typename Matrix::value_type, T>::value
                                     || boost::is_same<const typename Matrix::value_type,T>::value) );
            }

        template<typename Matrix2, typename U>
        explicit matrix_column_element_iterator(matrix_row_element_iterator<Matrix2,U> const& row_iter)
            : m(row_iter.m), row_pos(row_iter.row_pos), col_pos(row_iter.col_pos)
            {}
        
        template<typename Matrix2, typename U>
        matrix_column_element_iterator(matrix_column_element_iterator<Matrix2,U> const& r)
            : m(r.m), row_pos(r.row_pos), col_pos(r.col_pos)
            {}
    
    private:
        friend class boost::iterator_core_access;
        template <typename,typename> friend class matrix_column_element_iterator;

        value_type& dereference() const
        { return m->operator()(row_pos,col_pos); }

        // see comment for matrix_row_iterator::equal() and swap "row", "column"
        template <typename U>
        bool equal(matrix_column_element_iterator<Matrix,U> const& y) const
        {
            if(m == y.m && row_pos == y.row_pos)
                return true;
            else
                return false;
        }
        void increment()
        {
            ++(this->row_pos);
        }
        void decrement()
        {
            --(this->row_pos);
        }
        void advance(typename Matrix::difference_type n)
        {
            (this->row_pos) += n;
        }

        template <typename U>
        typename Matrix::difference_type distance_to(matrix_column_element_iterator<Matrix,U> const& z) const
        {
            return z.row_pos - row_pos;
        }
        
        typename Matrix::size_type row_pos;
        typename Matrix::size_type col_pos;
        Matrix* m;
};

#endif //BLAS_MATRIX_ITERATORS
