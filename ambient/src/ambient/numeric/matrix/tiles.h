#ifndef __AMBIENT_NUMERIC_MATRIX_H__
#define __AMBIENT_NUMERIC_MATRIX_H__

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

        explicit inline tiles();
        explicit inline tiles(size_type rows, size_type cols, value_type init_value = value_type()); 
        inline tiles(const tiles& m);
        tiles& operator = (const tiles& rhs); 
    public:
        inline size_type num_rows() const;
        inline size_type num_cols() const;
        inline scalar_type trace() const;
        inline void transpose();
        inline void conj();
        inline bool empty() const;          
        inline void swap(tiles& r);
        inline void resize(size_type rows, size_type cols); 
        inline void remove_rows(size_type i, size_type k = 1);
        inline void remove_cols(size_type j, size_type k = 1);
        inline Matrix& operator[] (size_type k);
        inline const Matrix& operator[] (size_type k) const;
        inline tiles& operator += (const tiles& rhs);
        inline tiles& operator -= (const tiles& rhs);
        template <typename T2> inline tiles& operator *= (const T2& t);
        template <typename T2> inline tiles& operator /= (const T2& t);
        inline value_type& operator() (size_type i, size_type j);
        inline const value_type& operator() (size_type i, size_type j) const;
        inline void load(alps::hdf5::archive & ar){};
        inline void save(alps::hdf5::archive & ar)const{};
    public:
        std::vector<Matrix*> data;
        size_t rows;
        size_t cols;
    };

} }

#endif
