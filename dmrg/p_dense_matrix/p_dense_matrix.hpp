#include "p_dense_matrix/p_dense_matrix.h"

namespace blas {

    template <typename T>
    p_dense_matrix<T>::p_dense_matrix(size_type rows = 0, size_type columns = 0, T init_value = T() )
    : rows(rows), cols(columns), lda(rows), sda(cols), action(NULL), proxy(false),
      data_scope(new T[rows*columns + ambient::get_bound(sizeof(T))]), 
      data(data_scope.get() + ambient::get_bound(sizeof(T)))
    {
        for(size_type i=0; i < rows*columns; i++) data[i] = init_value; // >_< // very slow //////////
//      memset(data, 0, rows*columns*sizeof(T)); // faster, though we don't need to initialize usually
    }

    template <typename T>
    p_dense_matrix<T>::p_dense_matrix(p_dense_matrix<T> const& m)
    : rows(m.rows), cols(m.cols), lda(m.lda), sda(m.sda), action(NULL), proxy(false),
      data_scope(new T[lda*sda + ambient::get_bound(sizeof(T))]),
      data(data_scope.get() + ambient::get_bound(sizeof(T)))
    {
        memcpy(this->data, m.data, this->lda*this->cols*sizeof(T));
    }

    template <typename T>
    p_dense_matrix<T>::~p_dense_matrix<T>(){
//        if(this->proxy) printf("Proxy matrix destroyed\n");
//        else printf("Non-proxy matrix destroyed\n");
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
                for(size_type i = this->cols; i-- >= 1;){
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
        for(size_type j = 0; j < this->cols; ++j)
            // for each column, copy the rows > i+k to k rows  up
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
    p_dense_matrix<T>& p_dense_matrix<T>::operator += (const p_dense_matrix& rhs){ return(*this = *this + rhs); }
    template <typename T>
    p_dense_matrix<T>& p_dense_matrix<T>::operator -= (const p_dense_matrix& rhs){ return (*this = *this - rhs); }
    template <typename T>
    template <typename T2>
    p_dense_matrix<T>& p_dense_matrix<T>::operator *= (const T2& t){ return (*this = *this * t); }
    template <typename T>
    template <typename T2>
    p_dense_matrix<T>& p_dense_matrix<T>::operator /= (const T2& t){ return (*this = *this * (1/t)); }


//////////////////////////////////// AMBIENT PART ////////////////////////////////////////////////////
    template <typename T>
    const ambient::p_profile p_dense_matrix<T>::profile() const { 
        if(proxy) return ambient::p_profile(this, "proxy"); 
        else return ambient::p_profile(this, "matrix"); 
    }
    template <typename T> // proxy object construction
    p_dense_matrix<T>::p_dense_matrix(const ambient::p_action* a): action(a), proxy(true), data_scope(new T()) { }
    template <typename T>
    p_dense_matrix<T>& p_dense_matrix<T>::operator = (p_dense_matrix const& rhs) // watch out of copying
    {
        this->action = new ambient::p_action('=', *this, rhs);
        return *this;
    }

    template <typename T>
    const p_dense_matrix<T> operator + (const p_dense_matrix<T>& a, const p_dense_matrix<T>& b){ return p_dense_matrix<T>(new ambient::p_action('+', a, b)); }
    template <typename T>
    const p_dense_matrix<T> operator - (const p_dense_matrix<T>& a, const p_dense_matrix<T>& b){ return p_dense_matrix<T>(new ambient::p_action('-', a, b)); }
    template<typename T>
    const p_dense_matrix<T> operator * (const p_dense_matrix<T>& lhs, const p_dense_matrix<T>& rhs){ return p_dense_matrix<T>(new ambient::p_action('*', lhs, rhs)); }
    template<typename T, typename T2>
    const p_dense_matrix<T> operator * (const p_dense_matrix<T>& m, const T2& t){ return p_dense_matrix<T>(new ambient::p_action('*', m, t)); }
    template<typename T, typename T2>
    const p_dense_matrix<T> operator * (const T2& t, const p_dense_matrix<T>& m){ return p_dense_matrix<T>(new ambient::p_action('*', t, m)); }
//////////////////////////////////// AMBIENT PART ////////////////////////////////////////////////////

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
