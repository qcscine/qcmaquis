#ifndef HDF5_STORAGE_H
#define HDF5_STORAGE_H

#include "utils/timings.h"

#include <boost/thread.hpp>

#include <deque>

//#define HDF5_STORAGE_VERBOSE

class Hdf5StorageMaster;

template<class Object>
class Hdf5Storage : public BaseStorage<Object>
{
public:
    Hdf5Storage(std::string const & object_path_, Hdf5StorageMaster * master_)
    : object_path(object_path_)
    , master(master_)
    {
    }
    
    boost::shared_ptr<BaseStorage<Object> > clone() const;
            
    void load(boost::shared_ptr<storage<Object> > ptr);
    void prefetch(boost::shared_ptr<storage<Object> > ptr);
    bool prefetch_barrier(boost::shared_ptr<storage<Object> > ptr);
    void store(boost::shared_ptr<storage<Object> > ptr);
    ~Hdf5Storage();
    
protected:
    std::string object_path;
    Hdf5StorageMaster * master;
    
    std::deque<boost::shared_ptr<boost::mutex> > prefetch_mutex;
    
    template<class Object_> friend class Hdf5ReadRequest_impl;
    template<class Object_> friend class Hdf5WriteRequest_impl;
};

class Hdf5Worker
{
public:
    Hdf5Worker(Hdf5StorageMaster * master_);
    
    void operator()();
    void execute();
private:
    Hdf5StorageMaster * master;
};

class Hdf5Request
{
public:
    virtual void operator()() = 0;
};

class Hdf5StorageMaster : public BaseStorageMaster
{
public:
    Hdf5StorageMaster(const char * fp)
    : file_path(fp)
    , last_id(0)
    , worker(this)
    , worker_thread(worker)
    , active(true)
    { }
    
    ~Hdf5StorageMaster()
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
    
    int flag() { return 1; }
    
    std::string get_path()
    {
        std::ostringstream oss;
        oss << "/" << last_id++;
        return oss.str();
    }
    
    template<class Object>
    boost::shared_ptr<BaseStorage<Object> > make()
    {
        return boost::shared_ptr<BaseStorage<Object> >(new Hdf5Storage<Object>(get_path(), this));
    }
    
    void notify()
    {
        worker_cond.notify_one();
    }
    
    void open_reading()
    {   
        oa.reset();
        
        if (ia.get() == NULL)
            ia.reset(new alps::hdf5::iarchive(file_path));
    }
    
    void open_writing()
    {   
        ia.reset();
        
        if (oa.get() == NULL)
            oa.reset(new alps::hdf5::oarchive(file_path));
    }
    
protected:
    std::string file_path;
    std::size_t last_id;
    
    boost::mutex h5mtx;
    
    std::deque<boost::shared_ptr<Hdf5Request> > requests;
    
    boost::condition_variable worker_cond;
//    boost::mutex worker_mut;
    
    boost::mutex deque_mutex;
    
    friend class Hdf5Worker;
    
    template<class Object>
    friend class Hdf5Storage;
    
    boost::shared_ptr<alps::hdf5::iarchive> ia;
    boost::shared_ptr<alps::hdf5::oarchive> oa;
    
    Hdf5Worker worker;
    boost::thread worker_thread;
    bool active;
    
    template<class Object_> friend class Hdf5ReadRequest_impl;
    template<class Object_> friend class Hdf5WriteRequest_impl;
    friend class Hdf5DeleteRequest;
};

template<class Object>
class Hdf5ReadRequest_impl : public Hdf5Request
{
public:
    Hdf5ReadRequest_impl(Hdf5StorageMaster * master_,
                         boost::shared_ptr<storage<Object> > ptr_,
                         boost::shared_ptr<boost::mutex> wait_mutex_)
    : master(master_)
    , ptr(ptr_)
    , wait_mutex(wait_mutex_)
    , wait_lock(*wait_mutex_)
    { }

    void operator()()
    {
        static Timer timer("hdf5_load");
        timer.begin();
        
        boost::lock_guard<boost::mutex> lock(master->h5mtx);
        
#ifdef HDF5_STORAGE_VERBOSE
        cerr << "Reading ";
#endif
  
        Hdf5Storage<Object> & store = dynamic_cast<Hdf5Storage<Object>&>(*ptr->storage_);
        
#ifdef HDF5_STORAGE_VERBOSE
        cerr << store.object_path << endl;
#endif
        
        try {
            master->open_reading();
            ptr->object_.reset(new Object());
            *master->ia >> alps::make_pvp(store.object_path, *ptr->object_);
        } catch (...) {
            cerr << "HDF5 threw an exception." << endl;
            exit(1);
        }
        
        timer.end();
    }
    
private:
    Hdf5StorageMaster * master;
    boost::shared_ptr<storage<Object> > ptr;
    boost::shared_ptr<boost::mutex> wait_mutex;
    boost::lock_guard<boost::mutex> wait_lock;
};

template<class Object>
class Hdf5WriteRequest_impl : public Hdf5Request
{
public:
    Hdf5WriteRequest_impl(Hdf5StorageMaster * master_,
                          boost::shared_ptr<storage<Object> > ptr_)
    : master(master_)
    , ptr(ptr_) { }
    
    void operator()()
    {
        static Timer timer("hdf5_store");
        timer.begin();
        
        boost::lock_guard<boost::mutex> lock(master->h5mtx);
        
#ifdef HDF5_STORAGE_VERBOSE
        cerr << "Writing ";
#endif
        
        Hdf5Storage<Object> & store = dynamic_cast<Hdf5Storage<Object>&>(*ptr->storage_);
        
#ifdef HDF5_STORAGE_VERBOSE
        cerr  << store.object_path << endl;
#endif
        
        try {
            master->open_writing();
            *master->oa << alps::make_pvp(store.object_path, *ptr->object_);
            ptr->object_.reset();
        } catch (...) {
            cerr << "HDF5 threw an exception." << endl;
            exit(1);
        }
        
#ifdef HDF5_STORAGE_VERBOSE
        cerr << "Done writing." << endl;
#endif
        
        timer.end();
    }
    
private:
    Hdf5StorageMaster * master;
    boost::shared_ptr<storage<Object> > ptr;
};

class Hdf5DeleteRequest : public Hdf5Request
{
public:
    Hdf5DeleteRequest(Hdf5StorageMaster * master_,
                      std::string path_)
    : master(master_)
    , object_path(path_) { }
    
    void operator()()
    {
        boost::lock_guard<boost::mutex> lock(master->h5mtx);
        
#ifdef HDF5_STORAGE_VERBOSE
        cerr << "Deleting " << object_path << endl;
#endif
        
        master->open_writing();
        if (master->oa->is_group(object_path))
            master->oa->delete_group(object_path);
    }
    
private:
    Hdf5StorageMaster * master;
    std::string object_path;
};

Hdf5Worker::Hdf5Worker(Hdf5StorageMaster * master_)
: master(master_)
{ }

void Hdf5Worker::operator()()
{
    while (true) {
        while (true) {
            boost::unique_lock<boost::mutex> lock(master->deque_mutex);
            if (!master->active)
                return;
            if (!master->requests.empty())
                break;
            master->worker_cond.wait(lock);
        }
        
        execute();
    }
}

void Hdf5Worker::execute()
{
    while (true) {
        boost::shared_ptr<Hdf5Request> req;
        {
            boost::lock_guard<boost::mutex> lock(master->deque_mutex);
            if (master->requests.empty())
                break;
            req = master->requests.front();
            master->requests.pop_front();
        }
        (*req)();
    }
}

template<> struct storage_master_type<1>
{
    typedef Hdf5StorageMaster type;
};

template<class Object>
boost::shared_ptr<BaseStorage<Object> > Hdf5Storage<Object>::clone() const
{
    return boost::shared_ptr<BaseStorage<Object> >(new Hdf5Storage<Object>(master->get_path(), master));
}

template<class Object>
void Hdf5Storage<Object>::load(boost::shared_ptr<storage<Object> > ptr)
{
#ifdef HDF5_STORAGE_VERBOSE
    cerr << "Load request " << object_path << endl;
#endif
    
    boost::shared_ptr<boost::mutex> wait_mutex(new boost::mutex());
    master->requests.push_back(boost::shared_ptr<Hdf5Request>(new Hdf5ReadRequest_impl<Object>(master, ptr, wait_mutex)));
    master->notify();
    
//    boost::lock_guard<boost::mutex> wait2(*wait_mutex);
    
    static Timer timer("load_wait");
    timer.begin();
    wait_mutex->lock();
    wait_mutex->unlock();
    timer.end();
}

template<class Object>
void Hdf5Storage<Object>::prefetch(boost::shared_ptr<storage<Object> > ptr)
{   
    assert(ptr->storage_.get() == this);
    
#ifdef HDF5_STORAGE_VERBOSE
    cerr << "Prefetch request " << object_path << endl;
#endif
    
    boost::shared_ptr<boost::mutex> wait_mutex(new boost::mutex());
    prefetch_mutex.push_back(wait_mutex);
    
    boost::lock_guard<boost::mutex> lock(master->deque_mutex);
    master->requests.push_back(boost::shared_ptr<Hdf5Request>(new Hdf5ReadRequest_impl<Object>(master, ptr, wait_mutex)));
    master->notify();
}

template<class Object>
void Hdf5Storage<Object>::store(boost::shared_ptr<storage<Object> > ptr)
{
    assert(ptr->storage_.get() == this);
    
#ifdef HDF5_STORAGE_VERBOSE
    cerr << "Store request " << object_path << endl;
#endif
    
    boost::lock_guard<boost::mutex> lock(master->deque_mutex);
    master->requests.push_back(boost::shared_ptr<Hdf5Request>(new Hdf5WriteRequest_impl<Object>(master, ptr)));
    master->notify();
}

template<class Object>
Hdf5Storage<Object>::~Hdf5Storage()
{
#ifdef HDF5_STORAGE_VERBOSE
    cerr << "Delete request " << object_path << endl;
#endif
    
    boost::lock_guard<boost::mutex> lock(master->deque_mutex);
    master->requests.push_back(boost::shared_ptr<Hdf5Request>(new Hdf5DeleteRequest(master, object_path)));
    master->notify();
}

template<class Object>
bool Hdf5Storage<Object>::prefetch_barrier(boost::shared_ptr<storage<Object> > ptr)
{
    assert(ptr->storage_.get() == this);
    
    if (prefetch_mutex.size() == 0)
        return false;
    
    static Timer timer("prefetch_barrier_wait");
    timer.begin();
    for (std::deque<boost::shared_ptr<boost::mutex> >::iterator it = prefetch_mutex.begin();
         it != prefetch_mutex.end(); ++it)
        boost::lock_guard<boost::mutex> wait(**it);
    prefetch_mutex.clear();
    timer.end();
    
    return true;
}
    
#endif
