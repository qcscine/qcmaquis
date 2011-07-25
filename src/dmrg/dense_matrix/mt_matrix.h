#ifndef MT_MATRIX_H
#define MT_MATRIX_H

#include <boost/thread.hpp>

#include "dense_matrix/dense_matrix.h"

#ifdef USE_GPU
#include "dense_matrix/gpu/m_matrix_gpu_functions.hpp"
#endif

#include "utils/data_collector.hpp"

#ifndef GEMM_MT_THRESHOLD
#define GEMM_MT_THRESHOLD 50
#endif
#ifndef GEMM_GPU_THRESHOLD
#define GEMM_GPU_THRESHOLD 100
#endif


template<typename T>
class mt_matrix;

class MtmRequest
{
public:
    virtual void operator()() = 0;
    virtual ~MtmRequest() { }
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
//        A_.wait(); B_.wait(); C_.wait();
        assert( A_.data_.num_cols() == B_.data_.num_rows() );
        static Timer timer("gemm requested");
        timer.begin();
#ifdef USE_GPU
        if (A_.data_.num_cols() > GEMM_GPU_THRESHOLD)
            gpu::matrix_matrix_multiply(A_.data_, B_.data_, C_.data_, 0.);
        else
            gemm(A_.data_, B_.data_, C_.data_);
#else
        gemm(A_.data_, B_.data_, C_.data_);
#endif
        timer.end();
    }
    
private:
    mt_matrix<T> const & A_;
    mt_matrix<T> const & B_;
    mt_matrix<T> & C_;
    boost::shared_ptr<boost::mutex> mutA_, mutB_, mutC_;
    boost::lock_guard<boost::mutex> lgA, lgB, lgC;
};

class MtmMaster;

class MtmWorker
{
public:
    MtmWorker(MtmMaster * m = NULL) : master(m) { };
    
    inline void operator()();
    inline void execute();
    
private:
    MtmMaster * master;
};

class MtmMaster
{
public:
    MtmMaster(int n = 1)
    {
        for (int i = 0; i < n; ++i)
            worker_threads.create_thread(MtmWorker(this));
    }
    
    ~MtmMaster ()
    {
        worker_threads.interrupt_all();
        notify();
        worker_threads.join_all();
    }
    
    void notify()
    {
        worker_cond.notify_one();
    }
    
    void push(boost::shared_ptr<MtmRequest> req)
    {
        boost::lock_guard<boost::mutex> lock(queue_mutex);
        queue.push_back(req);
        notify();
    }
    
private:
    std::list<boost::shared_ptr<MtmRequest> > queue;
    boost::mutex queue_mutex;
    
    boost::condition_variable worker_cond;
    
    template<typename T> friend class mt_matrix;
    friend class MtmWorker;
    
    boost::thread_group worker_threads;
};

#ifdef USE_MTM_MAIN
MtmMaster mtm_global_master;
#else
extern MtmMaster mtm_global_master;
#endif

void MtmWorker::operator()()
{
    while (true) {
        while (true) {
            boost::unique_lock<boost::mutex> lock(master->queue_mutex);
//            if (!master->active)
//                return;
            if (!master->queue.empty())
                break;
            master->worker_cond.wait(lock);
            boost::this_thread::interruption_point();
        }
        
        execute();
    }
}

void MtmWorker::execute()
{
    while (true) {
        boost::shared_ptr<MtmRequest> req;
        {
            boost::lock_guard<boost::mutex> lock(master->queue_mutex);
            if (master->queue.empty())
                break;
//            cerr << "Queue size: " << master->queue.size() << endl;
            req = master->queue.front();
            master->queue.pop_front();
        }
        (*req)();
    }
}

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
    
    explicit mt_matrix(size_type rows = 0, size_type cols = 0,
                       value_type val = value_type())
    : data_(rows, cols, val)
    , consistent(true)
    { }
    
    ~mt_matrix()
    {
        wait();
    }
    
    mt_matrix(mt_matrix const & rhs)
    {
        rhs.wait();
        data_ = rhs.data_;
        consistent = true;
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
    
    template<typename T2>
    mt_matrix& operator/=(T2 t)
    {
        wait();
        data_ /= t;
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
    
    void remove_cols(size_type i, difference_type k)
    {
        wait();
        data_.remove_cols(i, k);
    }
    
    void resize(size_type r, size_type c)
    {
        wait();
        data_.resize(r, c);
    }
    
#ifdef HAVE_ALPS_HDF5
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
#endif
    
    size_type num_rows() const { wait(); return data_.num_rows(); }
    size_type num_cols() const { wait(); return data_.num_cols(); }
    
    friend
    void svd(mt_matrix const & M,
             mt_matrix & U,
             mt_matrix & V,
             typename blas::diagonal_matrix<double> & S)
    {
        M.wait(); U.wait(); V.wait();
        svd(M.data_, U.data_, V.data_, S);
    }
    
    friend
    void gemm(mt_matrix const & A,
              mt_matrix const & B,
              mt_matrix & C)
    {
//        gemm(A.data_, B.data_, C.data_);
        
        A.wait(); B.wait(); C.wait();
        if (A.data_.num_cols() < GEMM_MT_THRESHOLD) {
            static Timer timer("single thread gemm");
            timer.begin();
            gemm(A.data_, B.data_, C.data_);
            timer.end();
        } else {
            boost::shared_ptr<boost::mutex> mta(new boost::mutex()), mtb(new boost::mutex()), mtc(new boost::mutex());
            {
                A.consistent = false;
                B.consistent = false;
                C.consistent = false;
                
                boost::lock_guard<boost::mutex> lock_listA(A.list_mutex), lock_listB(B.list_mutex), lock_listC(C.list_mutex);
                A.wait_mutices.push_back(mta);
                B.wait_mutices.push_back(mtb);
                C.wait_mutices.push_back(mtc);
            }
            
            DCOLLECTOR_ADD(gemm_collector, A.data_.num_cols())
            
            mtm_global_master.push(boost::shared_ptr<MtmRequest>(new MtmGemmRequest<T>(A, B, C,
                                                                                       mta, mtb, mtc)));
        }

//        cerr << "Waiting..." << endl;
//        A.wait(); B.wait(); C.wait();
//        cerr << "Done waiting." << endl;
    }
    
    friend
    void qr(mt_matrix const & M,
            mt_matrix & Q,
            mt_matrix & R)
    {
        M.wait(); Q.wait(); R.wait();
        qr(M.data_, Q.data_, R.data_);
    }
    
    friend
    void syev(mt_matrix const & M,
              mt_matrix & vecs,
              typename blas::diagonal_matrix<double> & S)
    {
        M.wait();
        vecs.wait();
        syev(M.data_, vecs.data_, S);
    }
    
    template<class Generator>
    void generate(Generator g)
    {
        wait();
        blas::generate(data_, g);
    }
    
    template<typename T2>
    friend mt_matrix operator*(mt_matrix m, T2 t)
    {
        m.wait();
        m *= t;
        return m;
    }
    
    template<typename T2>
    friend mt_matrix operator*(T2 t, mt_matrix m)
    {
        m.wait();
        m *= t;
        return m;
    }
    
private:
    slave_t data_;
    
    mutable boost::mutex list_mutex;
    mutable std::list<boost::shared_ptr<boost::mutex> > wait_mutices;
    mutable bool consistent;
    
    void wait() const
    {
        if (consistent)
            return;
        
        boost::lock_guard<boost::mutex> lock_list(list_mutex);
        for (std::list<boost::shared_ptr<boost::mutex> >::iterator it = wait_mutices.begin();
             it != wait_mutices.end(); ++it)
            boost::lock_guard<boost::mutex> lock(**it);
        wait_mutices.clear();
        
        consistent = true;
    }
    
    friend class MtmGemmRequest<T>;
};

template<typename T>
std::size_t num_rows(mt_matrix<T> const & m)
{
    return m.num_rows();
}

template<typename T>
std::size_t num_cols(mt_matrix<T> const & m)
{
    return m.num_cols();
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
    mt_matrix<T> ret(num_cols(m), num_rows(m));
    
    for (std::size_t r = 0; r < num_rows(m); ++r)
        for (std::size_t c = 0; c < num_cols(m); ++c)
            ret(c, r) = m(r, c);
    return ret;
}

template<typename T>
mt_matrix<T> conjugate(mt_matrix<T> m)
{
    for (std::size_t r = 0; r < num_rows(m); ++r)
        for (std::size_t c = 0; c < num_cols(m); ++c)
            m(r, c) = conj(m(r, c));
    return m;
}

template<typename T>
T trace(mt_matrix<T> const & m)
{
    T ret = T();
    for (std::size_t r = 0; r < std::min(num_rows(m), num_cols(m)); ++r)
        ret += m(r, r);
    return ret;
}

template<typename T>
std::ostream& operator<<(std::ostream & os, mt_matrix<T> const & m)
{
    for (std::size_t r = 0; r < num_rows(m); ++r) {
        for (std::size_t c = 0; c < num_cols(m); ++c)
            os << m(r, c);
        os << std::endl;
    }
    return os;
}

namespace blas
{
    template<typename T>
    struct associated_diagonal_matrix<mt_matrix<T> >
    {
        typedef typename associated_diagonal_matrix<dense_matrix<T> >::type type;
    };
    
    template<typename T>
    struct associated_vector<mt_matrix<T> >
    {
        typedef typename associated_vector<dense_matrix<T> >::type type;
    };
    
    template<typename T, class Generator>
    void generate(mt_matrix<T> & m, Generator g)
    {
        m.generate(g);
    }
}

#endif
