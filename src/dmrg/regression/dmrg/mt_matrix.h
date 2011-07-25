#ifndef MT_MATRIX_H
#define MT_MATRIX_H

#include <boost/thread.hpp>

#include "dense_matrix/dense_matrix.h"

class MtmRequest
{
public:
    virtual void operator()() = 0;
};

template<typename T>
class MtmGemmRequest : public MtmRequest
{
public:
    MtmGemmRequest(mt_matrix<T> const & A,
                   mt_matrix<T> const & B,
                   mt_matrix<T> & C,
                   boost::shared_ptr<boost::mutex> mutA,
                   boost::shared_ptr<boost::mutex> mutB,
                   boost::shared_ptr<boost::mutex> mutC)
    : A_(A), B_(B) , C_(C)
    , mutA_(mutA), mutB_(mutB), mutC_(mutC)
    , lgA(*mutA), lgB(*mutB), lgC(*mutC)
    { }
    
    void operator()()
    {
        gemm(A, B, C);
    }
    
private:
    mt_matrix<T> const & A, const & B, & C;
    boost::shared_ptr<boost::mutex> mutA_, mutB_, mutC_;
    boost::lock_guard<boost::mutex> lgA, lgB, lgC;
};

template<typename T>
class mt_matrix
{   
    typedef blas::dense_matrix<T> slave_t;
    
public:
    typedef typename slave_t::value_type value_type;
    typedef typename slave_t::reference reference;
    typedef typename slave_t::const_reference const_reference;
    typedef typename slave_t::size_type size_type;
    typedef typename slave_t::difference_type difference_type;
    
    typedef typename slave_t::element_iterator element_iterator;
    
    mt_matrix(size_type rows = 0, size_type columns = 0,
              value_type val = value_type())
    : data_(rows, columns, val)
    { }
    
    mt_matrix(mt_matrix const & rhs)
    {
        rhs.wait();
        data_ = rhs.data_;
    }
    
    void swap_with(mt_matrix & r)
    {
        wait(); r.wait();
        swap(this->data_, r.data_);
    }
    
    friend void swap(mt_matrix & x,
                     mt_matrix & y)
    {
        x.swap_with(y);
    }
    
    mt_matrix& operator=(mt_matrix r)
    {
        this->swap_with(r);
        return *this;
    }
    
    value_type & operator()(size_type i, size_type j)
    {
        wait();
        return data_(i, j);
    }
    
    value_type const & operator()(size_type i, size_type j) const
    {
        wait();
        return data_(i, j);
    }
    
    mt_matrix& operator+=(mt_matrix const & rhs)
    {
        wait(); rhs.wait();
        data_ += rhs.data_;
        return *this;
    }
    
    mt_matrix& operator-=(mt_matrix const & rhs)
    {
        wait(); rhs.wait();
        data_ -= rhs.data_;
        return *this;
    }
    
    template<typename T2>
    mt_matrix& operator*=(T2 t)
    {
        wait();
        data_ *= t;
        return *this;
    }
    
    std::pair<element_iterator,element_iterator> elements()
    {
        wait();
        return data_.elements();
    }
    
    void remove_rows(size_type i, difference_type k)
    {
        wait();
        data_.remove_rows(i, k);
    }
    
    void remove_columns(size_type i, difference_type k)
    {
        wait();
        data_.remove_columns(i, k);
    }
    
    void resize(size_type r, size_type c)
    {
        wait();
        data_.resize(r, c);
    }
    
    void load(alps::hdf5::archive & ar)
    {
        wait();
        data_.load(ar);
    }
    void save(alps::hdf5::archive & ar) const
    {
        wait();
        data_.save(ar);
    }
    
    size_type num_rows() const { return data_.num_rows(); }
    size_type num_columns() const { return data_.num_columns(); }
    
    friend
    void svd(mt_matrix const & M,
             mt_matrix & U,
             mt_matrix & V,
             typename blas::diagonal_matrix<double> & S)
    {
        svd(M.data_, U.data_, V.data_, S);
    }
    
    friend
    void gemm(mt_matrix const & A,
              mt_matrix const & B,
              mt_matrix & C)
    {
//        gemm(A.data_, B.data_, C.data_);
        boost::shared_ptr<boost::mutex> mta(new boost::mutex()), mtb(new boost::mutex()), mtc(new boost::mutex());
        A.wait_mutices.push_back(mta);
        B.wait_mutices.push_back(mtb);
        C.wait_mutices.push_back(mtc);
        
        queue.
    }
    
    friend
    void qr(mt_matrix const & M,
            mt_matrix & Q,
            mt_matrix & R)
    {
        qr(M.data_, Q.data_, R.data_);
    }
    
    friend
    void syev(mt_matrix const & M,
              mt_matrix & vecs,
              typename blas::diagonal_matrix<double> & S)
    {
        syev(M.data_, vecs.data_, S);
    }
    
private:
    slave_t data_;
    
    std::list<boost::shared_ptr<MtmRequest> > queue;
    boost::mutex queue_mutex;
    
    mutable std::list<boost::shared_ptr<boost::mutex> > wait_mutices;
    
    void wait() const
    {
        for (std::list<boost::shared_ptr<boost::mutex> >::iterator it = wait_mutices.begin();
             it != wait_mutices.end(); ++it)
            boost::lock_guard<boost::mutex> lock(**it);
        wait_mutices.clear();
    }
};

template<typename T>
std::size_t num_rows(mt_matrix<T> const & m)
{
    return m.num_rows();
}

template<typename T>
std::size_t num_columns(mt_matrix<T> const & m)
{
    return m.num_columns();
}

template<typename T>
std::pair<typename mt_matrix<T>::element_iterator, typename mt_matrix<T>::element_iterator>
elements(mt_matrix<T> & m)
{
    return m.elements();
}

template<typename T>
void resize(mt_matrix<T> & m,
            typename mt_matrix<T>::size_type r,
            typename mt_matrix<T>::size_type c)
{
    m.resize(r, c);
}

template<typename T>
mt_matrix<T> transpose(mt_matrix<T> const & m)
{
    mt_matrix<T> ret(num_columns(m), num_rows(m));
    
    for (std::size_t r = 0; r < num_rows(m); ++r)
        for (std::size_t c = 0; c < num_columns(m); ++c)
            ret(c, r) = m(r, c);
    return ret;
}

template<typename T>
mt_matrix<T> conjugate(mt_matrix<T> m)
{
    for (std::size_t r = 0; r < num_rows(m); ++r)
        for (std::size_t c = 0; c < num_columns(m); ++c)
            m(r, c) = conj(m(r, c));
    return m;
}

template<typename T>
T trace(mt_matrix<T> const & m)
{
    T ret = T();
    for (std::size_t r = 0; r < std::min(num_rows(m), num_columns(m)); ++r)
        ret += m(r, r);
    return ret;
}

template<typename T>
std::ostream& operator<<(std::ostream & os, mt_matrix<T> const & m)
{
    for (std::size_t r = 0; r < num_rows(m); ++r) {
        for (std::size_t c = 0; c < num_columns(m); ++c)
            os << m(r, c);
        os << endl;
    }
    return os;
}

#endif
