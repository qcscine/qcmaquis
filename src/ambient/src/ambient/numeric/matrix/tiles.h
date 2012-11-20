#ifndef __AMBIENT_NUMERIC_TILES_H__
#define __AMBIENT_NUMERIC_TILES_H__

namespace ambient { namespace numeric {

    template <class Matrix>
    class tiles {
    public:
        typedef typename Matrix::value_type  value_type;
        typedef typename Matrix::size_type   size_type;
        typedef typename Matrix::difference_type difference_type;
        typedef typename Matrix::real_type   real_type;
        typedef typename Matrix::scalar_type scalar_type;

        inline void* operator new (size_t);
        inline void operator delete (void* ptr);
        static inline tiles<Matrix> identity_matrix(size_type size);

        inline ~tiles();
        explicit inline tiles();
        explicit inline tiles(Matrix* a);
        explicit inline tiles(size_type rows, size_type cols, value_type init_value = value_type());
        inline tiles<subset_view<Matrix> > subset(size_type i, size_type j, size_type mt, size_type nt) const;
        inline tiles(const tiles& a);
        tiles& operator = (const tiles& rhs); 
        inline size_type num_rows() const;
        inline size_type num_cols() const;
        inline scalar_type trace() const;
        inline void transpose();
        inline void conj();
        inline bool empty() const;          
        inline void swap(tiles& r);
        inline void resize(size_type m, size_type n); 
        inline void remove_rows(size_type i, size_type k = 1);
        inline void remove_cols(size_type j, size_type k = 1);
        inline Matrix& tile(size_type i, size_type j);
        inline const Matrix& tile(size_type i, size_type j) const;
        inline Matrix& locate(size_type i, size_type j);
        inline const Matrix& locate(size_type i, size_type j) const;
        inline size_t addr(size_type i, size_type j) const;
        inline Matrix& operator[] (size_type k);
        inline const Matrix& operator[] (size_type k) const;
        template <class MatrixB> inline tiles& operator += (const tiles<MatrixB>& rhs);
        template <class MatrixB> inline tiles& operator -= (const tiles<MatrixB>& rhs);
        template <typename T2> inline tiles& operator *= (const T2& t);
        template <typename T2> inline tiles& operator /= (const T2& t);
        inline value_type& operator() (size_type i, size_type j);
        inline const value_type& operator() (size_type i, size_type j) const;
        inline void load(alps::hdf5::archive & ar){};
        inline void save(alps::hdf5::archive & ar)const{};
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
        inline Matrix& tile(size_type i, size_type j);
        inline const Matrix& tile(size_type i, size_type j) const;
        inline Matrix& operator[] (size_type k);
        inline const Matrix& operator[] (size_type k) const;
        inline tiles<subset_view<Matrix> > subset(size_type i, size_type j, size_type mt, size_type nt) const;
        template <class MatrixB> inline tiles& operator += (const tiles<MatrixB>& rhs);
        template <class MatrixB> inline tiles& operator -= (const tiles<MatrixB>& rhs);
        std::vector<subset_view<Matrix> > data;
        inline size_type num_rows() const;
        inline size_type num_cols() const;
        inline Matrix& locate(size_type i, size_type j);
        inline const Matrix& locate(size_type i, size_type j) const;
        inline size_t addr(size_type i, size_type j) const;
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

        inline void* operator new (size_t);
        inline void operator delete (void* ptr);

        explicit inline tiles(size_type size, value_type init_value = value_type()); 
        inline tiles(const tiles& a);
        tiles& operator = (const tiles& rhs); 
        inline ~tiles();
    public:
        inline size_type num_rows() const;
        inline size_type num_cols() const;
        inline void swap(tiles& r);
        inline void resize(size_type m, size_type n); 
        inline void remove_rows(size_type i, size_type k = 1);
        inline void remove_cols(size_type j, size_type k = 1);
        inline diagonal_matrix<T>& operator[] (size_type k);
        inline const diagonal_matrix<T>& operator[] (size_type k) const;
        inline value_type& operator() (size_type i, size_type j);
        inline const value_type& operator() (size_type i, size_type j) const;
    public:
        std::vector<diagonal_matrix<T>*> data;
        size_type size;
        size_type nt;
    };

} }

#endif
