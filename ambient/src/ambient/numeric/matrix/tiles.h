#ifndef AMBIENT_NUMERIC_TILES_H
#define AMBIENT_NUMERIC_TILES_H

namespace ambient { namespace numeric {

    template <class Matrix>
    class tiles {
    public:
        typedef typename Matrix::value_type  value_type;
        typedef typename Matrix::size_type   size_type;
        typedef typename Matrix::difference_type difference_type;
        typedef typename Matrix::real_type   real_type;
        typedef typename Matrix::scalar_type scalar_type;

        void* operator new (size_t);
        void operator delete (void* ptr);
        static tiles<Matrix> identity_matrix(size_type size);

       ~tiles();
        explicit tiles();
        explicit tiles(Matrix* a);
        explicit tiles(size_type rows, size_type cols, value_type init_value = value_type());
        tiles<subset_view<Matrix> > subset(size_type i, size_type j, size_type mt, size_type nt) const;
        tiles(const tiles& a);
        tiles& operator = (const tiles& rhs);
        size_type num_rows() const;
        size_type num_cols() const;
        scalar_type trace() const;
        void transpose();
        void conj();
        bool empty() const;          
        void swap(tiles& r);
        void resize(size_type m, size_type n); 
        void remove_rows(size_type i, size_type k = 1);
        void remove_cols(size_type j, size_type k = 1);
        Matrix& tile(size_type i, size_type j);
        const Matrix& tile(size_type i, size_type j) const;
        Matrix& locate(size_type i, size_type j);
        const Matrix& locate(size_type i, size_type j) const;
        size_t addr(size_type i, size_type j) const;
        Matrix& operator[] (size_type k);
        const Matrix& operator[] (size_type k) const;
        template <class MatrixB> operator tiles<MatrixB> () const;
        template <class MatrixB> tiles& operator += (const tiles<MatrixB>& rhs);
        template <class MatrixB> tiles& operator -= (const tiles<MatrixB>& rhs);
        template <typename T2> tiles& operator *= (const T2& t);
        template <typename T2> tiles& operator /= (const T2& t);
        value_type& operator() (size_type i, size_type j);
        const value_type& operator() (size_type i, size_type j) const;
        void load(alps::hdf5::archive & ar){};
        void save(alps::hdf5::archive & ar)const{};
    public:
        std::vector<Matrix*> data;
        size_type rows;
        size_type cols;
        size_type mt;
        size_type nt;
    };

    template <class Matrix>
    class tiles<subset_view<Matrix> >{
    public:
        typedef typename Matrix::size_type  size_type;
        typedef typename Matrix::value_type value_type;
        Matrix& tile(size_type i, size_type j);
        const Matrix& tile(size_type i, size_type j) const;
        Matrix& operator[] (size_type k);
        const Matrix& operator[] (size_type k) const;
        tiles<subset_view<Matrix> > subset(size_type i, size_type j, size_type mt, size_type nt) const;
        template <class MatrixB> tiles& operator += (const tiles<MatrixB>& rhs);
        template <class MatrixB> tiles& operator -= (const tiles<MatrixB>& rhs);
        std::vector<subset_view<Matrix> > data;
        size_type num_rows() const;
        size_type num_cols() const;
        Matrix& locate(size_type i, size_type j);
        const Matrix& locate(size_type i, size_type j) const;
        size_t addr(size_type i, size_type j) const;
        size_type rows;
        size_type cols;
        size_type mt;
        size_type nt;
    };

    template <typename T>
    class tiles<diagonal_matrix<T> > {
    public:
        typedef typename diagonal_matrix<T>::value_type  value_type;
        typedef typename diagonal_matrix<T>::size_type   size_type;
        typedef typename diagonal_matrix<T>::real_type   real_type;
        typedef typename diagonal_matrix<T>::scalar_type scalar_type;
        typedef typename diagonal_matrix<T>::difference_type difference_type;

        void* operator new (size_t);
        void operator delete (void* ptr);

        explicit tiles(size_type size, value_type init_value = value_type()); 
        tiles(const tiles& a);
        tiles& operator = (const tiles& rhs); 
       ~tiles();
    public:
        size_type num_rows() const;
        size_type num_cols() const;
        void swap(tiles& r);
        void resize(size_type m, size_type n); 
        void remove_rows(size_type i, size_type k = 1);
        void remove_cols(size_type j, size_type k = 1);
        diagonal_matrix<T>& operator[] (size_type k);
        const diagonal_matrix<T>& operator[] (size_type k) const;
        value_type& operator() (size_type i, size_type j);
        const value_type& operator() (size_type i, size_type j) const;
    public:
        std::vector<diagonal_matrix<T>*> data;
        size_type size;
        size_type nt;
    };

} }

#endif
