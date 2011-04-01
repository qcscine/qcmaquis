#ifndef DIAGONAL_MATRIX_H
#define DIAGONAL_MATRIX_H

namespace blas {

    template<class FullMatrixClass>
    struct associated_diagonal_matrix { };
    
	/**
	idea the container is a p_dense_matrix of one row and n columns,
	we need it to avoid dependency inside SVD kernels and other 
	*/	
    template<typename T>
    class p_diagonal_matrix
    {
    public:
        typedef T                       value_type;
        typedef T&                      reference;
        typedef T const&                const_reference;
        typedef std::size_t             size_type;
        typedef std::ptrdiff_t          difference_type;
        
        typedef typename std::vector<T>::iterator element_iterator;
        typedef typename std::vector<T>::const_iterator const_element_iterator;

		//to do
        template<class Vector>
        diagonal_matrix(Vector const & init)
        : data_(init.begin(), init.end()) { }

		//to do
        diagonal_matrix(std::size_t size = 0, T const & init = T())
        : data_(size, init) { }
        
        std::size_t num_rows() const { return data_.num_rows(); }
        std::size_t num_columns() const { return data_.num_rows(); }
        
        T const & operator[](std::size_t i) const { return data_(i,1); }
        T & operator[](std::size_t i) { return data_(i,1); }
        
        T const & operator()(std::size_t i, std::size_t j) const
        {
            assert(i == j);
            return data_(i,1);
        }
        T & operator()(std::size_t i, std::size_t j)
        {
            assert(i == j);
            return data_(i,1);
        }
        
        std::pair<element_iterator, element_iterator> elements()
        {
			int n = data.num_rows(); 
			return std::make_pair(data_.(1,1), data_.(n,1));
        }
        
        std::pair<const_element_iterator, const_element_iterator> elements() const
        {
 			int n = data.num_rows(); 
			return std::make_pair(data_.(1,1), data_.(n,1));
        }
        
        void remove_rows(std::size_t k, std::size_t n = 1)
        {
			remove_rows(k, k+n);

        }
        
        void remove_cols(std::size_t k, std::size_t n)
        {
			remove_rows(k, k+n);
		}
        
		//to do
        void resize(std::size_t r, std::size_t c, T v = T())
        {
            assert(r == c);
            data_.resize(r, v);
        }
        
    private:
		ambient::p_dense_matrix<T> data_;
    };
	
	//all free functions are inside p_diagonal
	// to do migrate all  class functions/constructors/destructor inside
	
  }

#endif
