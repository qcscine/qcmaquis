/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *               2011-2012 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef STREAM_STORAGE_H
#define STREAM_STORAGE_H

#include <boost/shared_ptr.hpp>
#include <boost/thread.hpp>

#include <iostream>
#include <fstream>
#include <deque>

#include "utils.hpp"
#include "utils/timings.h"


/********************
 * Forward declarations
 ********************/
template<class Matrix, class SymmGroup>
class block_matrix;
template<class Matrix, class SymmGroup>
class Boundary;

class StreamStorageMaster;


/********************
 * StreamStorage
 ********************/
class StreamStorage
{
public:
    StreamStorage(std::string const &, StreamStorageMaster*);
    
    StreamStorage(StreamStorage const & rhs);
    
    StreamStorage & operator=(StreamStorage rhs);
    
    enum status_t { Prefetching, Complete, Stored };
    status_t const& status() const { return status_; }
    status_t& status() { return status_; }
    
    StreamStorageMaster* master()
    { return master_; }
    
    // this operation does not have to be protected by a mutex, since ss.mutexes is only
    // modified by the calling thread
    void add_dependency(boost::shared_ptr<boost::mutex> const& new_mutex)
    { mutexes.push_back(new_mutex); }
    
    // wait for all dependencies to be finished
    void wait();
    
    ~StreamStorage();
    
protected:
    std::string object_path;
    StreamStorageMaster * master_;
    
    std::deque<boost::shared_ptr<boost::mutex> > mutexes;
    
    // friends
    template<class Object_> friend class StreamReadRequest_impl;
    template<class Object_> friend class StreamWriteRequest_impl;
    friend class StreamDeleteRequest;
    
    status_t status_;
};


/********************
 * StreamWorker
 ********************/
class StreamWorker
{
public:
    StreamWorker(StreamStorageMaster * master_);
    
    void operator()();
    void execute();
private:
    StreamStorageMaster * master;
};


/********************
 * StreamRequest
 ********************/
class StreamRequest
{
public:
    virtual void operator()() = 0;
};

template<class Object>
class StreamReadRequest_impl : public StreamRequest
{
public:
    StreamReadRequest_impl(StreamStorage * store_, Object * ptr_,
                           boost::shared_ptr<boost::mutex> wait_mutex_) { }
    
    void operator()() { }
};

template<class Object>
class StreamWriteRequest_impl : public StreamRequest
{
public:
    StreamWriteRequest_impl(StreamStorage *, Object *,
                            boost::shared_ptr<boost::mutex> wait_mutex_) { }
    
    void operator()() { }
};


/********************
 * StreamStorageMaster
 ********************/
class StreamStorageMaster
{
public:
    typedef StreamStorage Storage;
    
    StreamStorageMaster(std::string fp, bool enable=false);
    ~StreamStorageMaster();
    StreamStorage child();
    
    std::string get_path() const;
    void print_size() const;
    
    void add_req(boost::shared_ptr<StreamRequest> const& req);
    void notify();
    void sync();
    
    std::string & get_base_path() {return base_path;}
protected:
    std::string base_path;
    mutable std::size_t last_id;
    
    boost::mutex deque_mutex;
    std::deque<boost::shared_ptr<StreamRequest> > requests;
    boost::condition_variable worker_cond;
    
    bool active;
    StreamWorker worker;
    boost::thread worker_thread;
    
    friend class StreamWorker;
    friend class StreamStorage;
    
    template<class Object_> friend class StreamReadRequest_impl;
    template<class Object_> friend class StreamWriteRequest_impl;
    friend class StreamDeleteRequest;
};


/********************
 * Stream__Request_impl for Boundary<Matrix, SymmGroup>
 ********************/

template<class Matrix, class SymmGroup>
class StreamReadRequest_impl<Boundary<Matrix, SymmGroup> >
: public StreamRequest
{
    typedef Boundary<Matrix, SymmGroup> Object;
    
public:
    StreamReadRequest_impl(StreamStorage * store_,
                           Object * ptr_,
                           boost::shared_ptr<boost::mutex> wait_mutex_)
    : store(store_)
    , ptr(ptr_)
    , wait_mutex(wait_mutex_)
    , wait_lock(*wait_mutex_)
    { }
    
    void operator()()
    {
        std::string fp = store->master()->base_path + store->object_path;
        std::ifstream ifs(fp.c_str(), std::ifstream::binary);
        if (!ifs) {
            maquis::cerr << "File not found in StreamReadRequest!" << std::endl;
            exit(1);
        }
        
        Object & o = *ptr;
        
        for (std::size_t i = 0; i < o.aux_dim(); ++i)
            for (std::size_t k = 0; k < o[i].n_blocks(); ++k) {
                o[i][k] = Matrix(o[i].left_basis()[k].second,
                                 o[i].right_basis()[k].second);
                ifs.read((char*)(&o[i][k](0,0)),
                         o[i].left_basis()[k].second*
                         o[i].right_basis()[k].second*
                         sizeof(typename Matrix::value_type)/sizeof(char));
            }
        
        ifs.close();
    }
    
private:
    StreamStorage * store;
    Object * ptr;
    boost::shared_ptr<boost::mutex> wait_mutex;
    boost::lock_guard<boost::mutex> wait_lock;
};

template<class Matrix, class SymmGroup>
class StreamWriteRequest_impl<Boundary<Matrix, SymmGroup> >
: public StreamRequest
{
    typedef Boundary<Matrix, SymmGroup> Object;
    
public:
    StreamWriteRequest_impl(StreamStorage * store_,
                            Object * ptr_,
                            boost::shared_ptr<boost::mutex> wait_mutex_)
    : store(store_)
    , ptr(ptr_)
    , wait_mutex(wait_mutex_)
    , wait_lock(*wait_mutex_)
    { }
    
    void operator()()
    {
        std::string fp = store->master()->base_path + store->object_path;
        std::ofstream of(fp.c_str(), std::ofstream::binary);
        
        Object & o = *ptr;
        
        for (std::size_t i = 0; i < o.aux_dim(); ++i)
            for (std::size_t k = 0; k < o[i].n_blocks(); ++k) {
                // workaround until capacity business in alps::numeric::matrix is sorted out
                for (std::size_t c = 0; c < num_cols(o[i][k]); ++c)
                    of.write((char*)(&o[i][k](0, c)),
                             num_rows(o[i][k]) *
                             sizeof(typename Matrix::value_type)/sizeof(char));
                o[i][k] = Matrix();
            }
        
        of.close();
    }
    
private:
    StreamStorage * store;
    Object * ptr;
    boost::shared_ptr<boost::mutex> wait_mutex;
    boost::lock_guard<boost::mutex> wait_lock;
};


/********************
 * StreamWriteRequest_impl for block_matrix<Matrix, SymmGroup>
 ********************/

template<class Matrix, class SymmGroup>
class StreamReadRequest_impl<block_matrix<Matrix, SymmGroup> >
: public StreamRequest
{
    typedef block_matrix<Matrix, SymmGroup> Object;
    
public:
    StreamReadRequest_impl(StreamStorage * store_,
                           Object * ptr_,
                           boost::shared_ptr<boost::mutex> wait_mutex_)
    : store(store_)
    , ptr(ptr_)
    , wait_mutex(wait_mutex_)
    , wait_lock(*wait_mutex_)
    { }
    
    void operator()()
    {
        std::string fp = store->master()->base_path + store->object_path;
        std::ifstream ifs(fp.c_str(), std::ifstream::binary);
        if (!ifs) {
            maquis::cerr << "File not found in StreamReadRequest!" << std::endl;
            exit(1);
        }
        
        Object & o = *ptr;
        
        for (std::size_t k = 0; k < o.n_blocks(); ++k) {
            o[k] = Matrix(o.left_basis()[k].second,
                          o.right_basis()[k].second);
            ifs.read((char*)(&o[k](0,0)),
                     o.left_basis()[k].second*
                     o.right_basis()[k].second*
                     sizeof(typename Matrix::value_type)/sizeof(char));
        }
        
        ifs.close();
    }
    
private:
    StreamStorage * store;
    Object * ptr;
    boost::shared_ptr<boost::mutex> wait_mutex;
    boost::lock_guard<boost::mutex> wait_lock;
};


template<class Matrix, class SymmGroup>
class StreamWriteRequest_impl<block_matrix<Matrix, SymmGroup> >
: public StreamRequest
{
    typedef block_matrix<Matrix, SymmGroup> Object;
    
public:
    StreamWriteRequest_impl(StreamStorage * store_,
                            Object * ptr_,
                            boost::shared_ptr<boost::mutex> wait_mutex_)
    : store(store_)
    , ptr(ptr_)
    , wait_mutex(wait_mutex_)
    , wait_lock(*wait_mutex_)
    { }
    
    void operator()()
    {
        std::string fp = store->master()->base_path + store->object_path;
        std::ofstream of(fp.c_str(), std::ofstream::binary);
        
        Object & o = *ptr;
        
        for (std::size_t k = 0; k < o.n_blocks(); ++k) {
            // workaround until capacity business in alps::numeric::matrix is sorted out
            for (std::size_t c = 0; c < num_cols(o[k]); ++c)
                of.write((char*)(&o[k](0, c)),
                         num_rows(o[k]) *
                         sizeof(typename Matrix::value_type)/sizeof(char));
            o[k] = Matrix();
        }
        
        of.close();
    }
    
private:
    StreamStorage * store;
    Object * ptr;
    boost::shared_ptr<boost::mutex> wait_mutex;
    boost::lock_guard<boost::mutex> wait_lock;
};

class StreamDeleteRequest : public StreamRequest
{
public:
    StreamDeleteRequest(StreamStorage * store)
    : file_path(store->master()->base_path + store->object_path) { }
    
    void operator()()
    {
        remove(file_path.c_str());
    }
    
private:
    std::string file_path;
};


/********************
 * storage
 ********************/
namespace storage {
    template<class T>
    void prefetch(T & o, StreamStorage & ss)
    {
        if (ss.master()->get_base_path().size() == 0)
            return;
        
        if (ss.status() == StreamStorage::Prefetching)
            return;
        if (ss.status() == StreamStorage::Complete)
            return;
        
        ss.status() = StreamStorage::Prefetching;
        
        boost::shared_ptr<boost::mutex> dependency_mutex(new boost::mutex());
        boost::shared_ptr<StreamRequest> req(new StreamReadRequest_impl<T>(&ss, &o, dependency_mutex));
        
        ss.master()->add_req(req);
        ss.add_dependency(dependency_mutex);
        ss.master()->notify();
    }
    
    template<class T>
    void load(T & o, StreamStorage & ss)
    {
        if (ss.master()->get_base_path().size() == 0)
            return;
        
        if (ss.status() == StreamStorage::Complete)
            return;
        else if (ss.status() == StreamStorage::Prefetching)
            ss.wait();
        else if (ss.status() == StreamStorage::Stored) {
            prefetch(o, ss);
            ss.wait();
        } else
            throw std::runtime_error("WTF?");
        assert(ss.status() == StreamStorage::Complete);
    }
    
    template<class T>
    void store(T & o, StreamStorage & ss)
    {
        if (ss.master()->get_base_path().size() == 0)
            return;
        
        if (ss.status() == StreamStorage::Stored)
            return;
        
        boost::shared_ptr<boost::mutex> dependency_mutex(new boost::mutex());
        boost::shared_ptr<StreamRequest> req(new StreamWriteRequest_impl<T>(&ss, &o, dependency_mutex));
        
        ss.master()->add_req(req);
        ss.add_dependency(dependency_mutex);
        ss.master()->notify();
        
        ss.status() = StreamStorage::Stored;
    }
    
    inline void reset(StreamStorage & ss)
    {
        if (ss.master()->get_base_path().size() == 0)
            return;
        
        ss.wait();
        assert(ss.status() == StreamStorage::Complete);
    }
    
};

#ifdef USE_AMBIENT
#include "stream_storage_ambient.h" // the specialization for ****_impl function object
#endif

            
#endif
