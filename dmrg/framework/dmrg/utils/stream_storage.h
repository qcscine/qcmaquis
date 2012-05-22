/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
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

template<class Matrix, class SymmGroup>
class block_matrix; 
template<class Matrix, class SymmGroup>
class Boundary; 


class NoopStorage { };
  
class NoopStorageMaster {
 public:
    typedef NoopStorage Storage;
    NoopStorage child() { return NoopStorage(); }
    void sync() const { }
  
 };

class StreamStorageMaster;

class StreamStorage
{
public:
    StreamStorage(std::string const & object_path_,
                  StreamStorageMaster * master_);
    
    StreamStorage(StreamStorage const & rhs);
    
    StreamStorage & operator=(StreamStorage rhs)
    {
        std::swap(*this, rhs);
        return *this;
    }
    
    friend struct storage;
    
    enum status_t { Prefetching, Complete, Stored };
    status_t status() { return status_; }
    
    ~StreamStorage();
    
protected:
    std::string object_path;
    StreamStorageMaster * master;
    
    std::deque<boost::shared_ptr<boost::mutex> > mutexes;
    
    template<class Object_> friend class StreamReadRequest_impl;
    template<class Object_> friend class StreamWriteRequest_impl;
    friend class StreamDeleteRequest;
    
    void wait()
    {
        for (std::deque<boost::shared_ptr<boost::mutex> >::iterator it = mutexes.begin();
             it != mutexes.end(); ++it)
            boost::lock_guard<boost::mutex> lock(**it);
        mutexes.clear();
        status_ = Complete;
    }
    
    status_t status_;
};

class StreamWorker
{
public:
    StreamWorker(StreamStorageMaster * master_);
    
    void operator()();
    void execute();
private:
    StreamStorageMaster * master;
};

class StreamRequest
{
public:
    virtual void operator()() = 0;
};

class StreamStorageMaster
{
public:
    typedef StreamStorage Storage;
    
    StreamStorageMaster(std::string fp, bool enable=false);
    
    ~StreamStorageMaster()
    {
        {
            boost::lock_guard<boost::mutex> lock(deque_mutex);
            active = false;
        }
        active = false;
        notify();
        worker_thread.join();
        worker.execute();
    }
    
    std::string get_path()
    {
        std::ostringstream oss;
        oss << "/" << last_id++;
        return oss.str();
    }
    
    void print_size()
    {
        if (base_path.size() == 0)
            return;
        
        maquis::cout << "Storage directory size: "; std::cout.flush();
        std::ostringstream oss;
        oss << "du -skh " << base_path << " | awk '{print $1}'";
        system(oss.str().c_str());
    }
    
    void sync()
    {
        {
            boost::lock_guard<boost::mutex> lock(deque_mutex);
            active = false;
        }
        notify();
        worker_thread.join();
        worker.execute();
        {
            boost::lock_guard<boost::mutex> lock(deque_mutex);
            active = true;
        }
        worker_thread = boost::thread(worker);
    }
    
    StreamStorage child()
    {
        return StreamStorage(get_path(), this);
    }
    
    void notify()
    {
        worker_cond.notify_one();
    }
    
protected:
    std::string base_path;
    std::size_t last_id;
    
    std::deque<boost::shared_ptr<StreamRequest> > requests;
    
    boost::condition_variable worker_cond;
    
    boost::mutex deque_mutex;
    
    friend class StreamWorker;
    friend class StreamStorage;
    
    bool active;
    StreamWorker worker;
    boost::thread worker_thread;
    
    template<class Object_> friend class StreamReadRequest_impl;
    template<class Object_> friend class StreamWriteRequest_impl;
    friend class StreamDeleteRequest;
    
    friend struct storage;
};

template<class Object>
class StreamReadRequest_impl : public StreamRequest
{
public:
    StreamReadRequest_impl(StreamStorage * store_,
                           Object * ptr_,
                           boost::shared_ptr<boost::mutex> wait_mutex_) { }
    
    void operator()() { }
};

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
        if (store->master->base_path.size() == 0)
            return;
        
        std::string fp = store->master->base_path + store->object_path;
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
        if (store->master->base_path.size() == 0)
            return;
        
        std::string fp = store->master->base_path + store->object_path;
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

template<class Object>
class StreamWriteRequest_impl : public StreamRequest
{
public:
    StreamWriteRequest_impl(StreamStorage *,
                            Object *) { }
    
    void operator()()
    { }
};

template<class Matrix, class SymmGroup>
class StreamWriteRequest_impl<Boundary<Matrix, SymmGroup> >
: public StreamRequest
{
    typedef Boundary<Matrix, SymmGroup> Object;
    
public:
    StreamWriteRequest_impl(StreamStorage * store_,
                            Object * ptr_)
    : store(store_)
    , ptr(ptr_) { }
    
    void operator()()
    {
        if (store->master->base_path.size() == 0)
            return;
        
        std::string fp = store->master->base_path + store->object_path;
        std::ofstream of(fp.c_str(), std::ofstream::binary);
        
        Object & o = *ptr;
        
        for (std::size_t i = 0; i < o.aux_dim(); ++i)
            for (std::size_t k = 0; k < o[i].n_blocks(); ++k) {
                // workaround until capacity business in dense_matrix is sorted out
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
};

template<class Matrix, class SymmGroup>
class StreamWriteRequest_impl<block_matrix<Matrix, SymmGroup> >
: public StreamRequest
{
    typedef block_matrix<Matrix, SymmGroup> Object;
    
public:
    StreamWriteRequest_impl(StreamStorage * store_,
                            Object * ptr_)
    : store(store_)
    , ptr(ptr_) { }
    
    void operator()()
    {
        
        if (store->master->base_path.size() == 0)
            return;
        
        std::string fp = store->master->base_path + store->object_path;
        std::ofstream of(fp.c_str(), std::ofstream::binary);
        
        Object & o = *ptr;
        
        for (std::size_t k = 0; k < o.n_blocks(); ++k) {
            // workaround until capacity business in dense_matrix is sorted out
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
};

class StreamDeleteRequest : public StreamRequest
{
public:
    StreamDeleteRequest(StreamStorage * store)
    : file_path(store->master->base_path + store->object_path) { }
    
    void operator()()
    {
        remove(file_path.c_str());
    }
    
private:
    std::string file_path;
};


struct storage {
    template<class T>
    static void prefetch(T & o,
                         StreamStorage & ss)
    {
        if (ss.status() == StreamStorage::Prefetching)
            return;
        if (ss.status() == StreamStorage::Complete)
            return;
        
        ss.status_ = StreamStorage::Prefetching;
        
        boost::shared_ptr<boost::mutex> new_mutex(new boost::mutex());
        boost::shared_ptr<StreamRequest> req(new StreamReadRequest_impl<T>(&ss, &o, new_mutex));
        
        {
            boost::lock_guard<boost::mutex> lock(ss.master->deque_mutex);
            ss.master->requests.push_back(req);
        }
        
        // this operation does not have to be protected by a mutex, since ss.mutexes is only
        // modified by the calling thread
        ss.mutexes.push_back(new_mutex);
        ss.master->notify();
    }
    
    template<class T>
    static void load(T & o,
                     StreamStorage & ss)
    {
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
    static void store(T & o,
                      StreamStorage & ss)
    {
        if (ss.status() == StreamStorage::Stored)
            return;
        
        boost::shared_ptr<StreamRequest> req(new StreamWriteRequest_impl<T>(&ss, &o));
        
        {
            boost::lock_guard<boost::mutex> lock(ss.master->deque_mutex);
            ss.master->requests.push_back(req);
        }
        ss.master->notify();
        
        ss.status_ = StreamStorage::Stored;
    }
    
    static void reset(StreamStorage & ss)
    {
        ss.wait();
        ss.status_ = StreamStorage::Complete;
    }

/* clean me up! */
 template<class T> static void store(T &, NoopStorage &) { }
 template<class T> static void prefetch(T &, NoopStorage &) { }
 template<class T> static void load(T &, NoopStorage &) { }
 static void reset(NoopStorage &) { }
};


#endif
