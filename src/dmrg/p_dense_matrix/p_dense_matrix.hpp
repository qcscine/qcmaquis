#include "p_dense_matrix/p_dense_matrix.h"

namespace blas {

    template <typename T>
    p_dense_matrix<T>::p_dense_matrix(size_type rows = 0, size_type columns = 0, T init_value = T() )
    : rows(rows), cols(columns), lda(rows), sda(cols), 
      data_scope(new T[rows*columns + ambient::get_bound(sizeof(T))]), 
      data(data_scope.get() + ambient::get_bound(sizeof(T)))
    {
        memset(data, 0, rows*columns*sizeof(T));
    }

    template <typename T>
    p_dense_matrix<T>::p_dense_matrix(p_dense_matrix<T> const& m)
    : rows(m.rows), cols(m.cols), lda(m.lda), sda(m.sda),
      data_scope(new T[lda*sda + ambient::get_bound(sizeof(T))]),
      data(data_scope.get() + ambient::get_bound(sizeof(T)))
    {
        memcpy(this->data, m.data, this->lda*this->cols*sizeof(T));
    }

    template <typename T>
    void p_dense_matrix<T>::swap(p_dense_matrix & r)
    {
        this->data_scope.swap(r.data_scope);
        std::swap(this->data, r.data);
        std::swap(this->rows, r.rows);
        std::swap(this->cols, r.cols);
        std::swap(this->lda,r.lda);
        std::swap(this->sda,r.sda);
    }

#ifdef HAVE_ALPS_HDF5
    template <typename T>
    inline void p_dense_matrix<T>::serialize(alps::hdf5::iarchive & ar)
    {
        ar >> alps::make_pvp("rows", rows);
        ar >> alps::make_pvp("cols", cols);
        ar >> alps::make_pvp("lda", lda);
        ar >> alps::make_pvp("values", data);
    }
    template <typename T>
    inline void p_dense_matrix<T>::serialize(alps::hdf5::oarchive & ar) const
    {
        ar << alps::make_pvp("rows", rows);
        ar << alps::make_pvp("cols", cols);
        ar << alps::make_pvp("lda", lda);
        ar << alps::make_pvp("values", data);
    }
#endif	

    template <typename T>
    inline const bool p_dense_matrix<T>::empty() const { return (this->rows == 0 || this->cols == 0); }

    template <typename T>
    inline const std::size_t p_dense_matrix<T>::num_rows() const { return this->rows; }

    template <typename T>
    inline const std::size_t p_dense_matrix<T>::num_columns() const { return this->cols; }

    template <typename T>
    inline const std::ptrdiff_t p_dense_matrix<T>::stride1() const { return 1; }

    template <typename T>
    inline const std::ptrdiff_t p_dense_matrix<T>::stride2() const { return this->lda; }

    template <typename T>
    void p_dense_matrix<T>::clear(){ this->rows = this->cols = 0; }

    template <typename T>
    void p_dense_matrix<T>::reserve(size_type rows, size_type cols)
    {
        if(rows > this->lda || cols > this->sda){
                size_type lda = std::max(rows,this->lda)*3/2; // reservation
                size_type sda = std::max(cols,this->sda)*3/2; // reservation
                boost::scoped_ptr<T> tmp(new T[ambient::get_bound(sizeof(T))+lda*sda]);
                this->data_scope.swap(tmp);
                this->data = this->data_scope.get() + ambient::get_bound(sizeof(T));
                memcpy(this->data_scope.get(), tmp.get(), 
                       sizeof(T)*(ambient::get_bound(sizeof(T))+this->lda*this->cols)); // restoring original data
                for(size_type i = this->cols-1; i >= 0; i--){
                    memmove(this->data+lda*i, this->data+this->lda*i, sizeof(T)*this->rows);
                }
                this->lda = lda;
                this->sda = sda;
        }
    }

    template <typename T>
    void p_dense_matrix<T>::resize(size_type rows, size_type cols)
    {
        reserve(rows,cols);
        if(rows > this->rows){
            for(size_type i = 0; i < this->cols; i++){
                memset(this->data+this->lda*i+this->rows, 0, sizeof(T)*(rows - this->rows));
            }
        }
        if(cols > this->cols){
            for(size_type i = this->cols; i < cols; i++){
                memset(this->data+this->lda*i, 0, sizeof(T)*rows);
            }
        }
        this->rows=rows;
        this->cols=cols;
    }

    template <typename T>
    void p_dense_matrix<T>::remove_rows(size_type i, difference_type k = 1)
    {
        assert( i+k <= this->rows );
        // for each column, copy the rows > i+k to k rows  up
        for(size_type j = 0; j < this->cols; ++j)
            memmove(&this->data[this->lda*j + i], &this->data[this->lda*j + i + k], sizeof(T)*(this->rows-i-k));
        this->rows -= k;
    }

    template <typename T>
    void p_dense_matrix<T>::remove_columns(size_type j, difference_type k = 1)
    {
        assert( j+k <= this->cols );
        this->cols -= k;
    }

    template <typename T>
    void p_dense_matrix<T>::inplace_conjugate()
    {
        assert(false);
        // TODO
    }

    template <typename T>
    inline T& p_dense_matrix<T>::operator()(const size_type i, const size_type j)
    {
        assert(i < this->rows);
        assert(j < this->cols);
        return this->data[i+j*this->lda];
    }

    template <typename T>
    inline T const& p_dense_matrix<T>::operator()(const size_type i, const size_type j) const 
    {
        assert((i < this->rows) && (j < this->cols));
        return this->data[i+j*this->lda];
    }

    template <typename T>
    p_dense_matrix<T>& p_dense_matrix<T>::operator = (p_dense_matrix<T> rhs) // watch out of copying
    {
        this->swap(rhs);
        return *this;
    }
    template <typename T>
    template <typename L, typename R>
    p_dense_matrix<T>& p_dense_matrix<T>::operator = (const p_action<L,R> & a) 
    {
// parse our p_actions here !!!

        return *this;
    }

    template <typename T>
    p_dense_matrix<T>& p_dense_matrix<T>::operator += (p_dense_matrix const& rhs){ return(*this = *this + rhs); }
    template <typename T>
    template <typename L, typename R>
    p_dense_matrix<T>& p_dense_matrix<T>::operator += (p_action<L,R> const& a){ return (*this = *this + a); }
    template <typename T>
    p_dense_matrix<T>& p_dense_matrix<T>::operator -= (p_dense_matrix const& rhs){ return (*this = *this - rhs); }
    template <typename T>
    template <typename L, typename R>
    p_dense_matrix<T>& p_dense_matrix<T>::operator -= (p_action<L,R> const& a){ return (*this = *this - a); }
    template <typename T>
    template <typename T2>
    p_dense_matrix<T>& p_dense_matrix<T>::operator *= (T2 const& t){ return (*this = *this * t); }
    template <typename T>
    template <typename L, typename R>
    p_dense_matrix<T>& p_dense_matrix<T>::operator *= (p_action<L,R> const& a){ return (*this = *this * a); }
    template <typename T>
    template <typename T2>
    p_dense_matrix<T>& p_dense_matrix<T>::operator /= (T2 const& t){ return *this *= (1/t); }
}

namespace blas {
    template <typename T>
    const p_action< p_dense_matrix<T>, p_dense_matrix<T> > operator + (p_dense_matrix<T> const& a, p_dense_matrix<T> const& b){ return p_action< p_dense_matrix<T>, p_dense_matrix<T> >('+', a, b); }
    template <typename T, typename L, typename R>
    const p_action< p_dense_matrix<T>, p_action<L, R> > operator + (p_dense_matrix<T> const& a, p_action<L,R> const& b){ return p_action< p_dense_matrix<T>,p_action<L,R> >('+', a, b); } 
    template <typename T, typename L, typename R>
    const p_action< p_action<L, R>, p_dense_matrix<T> > operator + (p_action<L,R> const& a, p_dense_matrix<T> const& b){ return p_action< p_action<L,R>,p_dense_matrix<T> >('+', a, b); }
    template <typename T>
    const p_action< p_dense_matrix<T>, p_dense_matrix<T> > operator - (p_dense_matrix<T> const& a, p_dense_matrix<T> const& b){ return p_action< p_dense_matrix<T>,p_dense_matrix<T> >('-', a, b); }
    template <typename T, typename L, typename R>
    const p_action< p_dense_matrix<T>, p_action<L, R> > operator - (p_dense_matrix<T> const& a, p_action<L,R> const& b){ return p_action< p_dense_matrix<T>,p_action<L,R> >('-', a, b); } 
    template <typename T, typename L, typename R>                                                                                                                                          
    const p_action< p_action<L, R>, p_dense_matrix<T> > operator - (p_action<L,R> const& a, p_dense_matrix<T> const& b){ return p_action< p_action<L,R>,p_dense_matrix<T> >('-', a, b); }
    template<typename T>
    const p_action< p_dense_matrix<T>, p_dense_matrix<T> > operator * (p_dense_matrix<T> const& lhs, p_dense_matrix<T> const& rhs){ return p_action< p_dense_matrix<T>,p_dense_matrix<T> >('*', lhs, rhs); }
    template<typename T, typename T2>
    const p_action< p_dense_matrix<T>, T2 > operator * (p_dense_matrix<T> const& m, T2 const& t){ return p_action< p_dense_matrix<T>, T2 >('*', m, t); }
    template<typename T, typename T2>
    const p_action< T2, p_dense_matrix<T> >operator * (T2 const& t, p_dense_matrix<T> const& m){ return p_action< T2, p_dense_matrix<T> >('*', t, m); }
    template <typename T, typename T2, typename L, typename R>
    const p_action< T2, p_action<L, R> > operator * (T2 const& t, p_action<L,R> const& a){ return p_action< T2, p_action<L,R> >('*', t, a); }
    template <typename T, typename T2, typename L, typename R>
    const p_action< p_action<L, R>, T2 > operator * (p_action<L,R> const& a, T2 const& t){ return p_action< p_action<L,R>, T2 >('*', a, t); }
    template <typename T>
    std::ostream& operator << (std::ostream& o, p_dense_matrix<T> const& m)
    {
        for(typename p_dense_matrix<T>::size_type i=0; i< m.num_rows(); ++i)
        {
            for(typename p_dense_matrix<T>::size_type j=0; j < m.num_columns(); ++j)
                o<<m(i,j)<<" ";
            o<<std::endl;
        }
        return o;
    }
}
