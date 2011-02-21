#ifndef HDF5_STORAGE_H
#define HDF5_STORAGE_H

#include <utils/timings.h>

#include <boost/thread.hpp>

class Hdf5StorageMaster;

template<class Object>
class Hdf5Storage : public BaseStorage<Object>
{
public:
    Hdf5Storage(std::string const & object_path_, Hdf5StorageMaster * master_)
    : object_path(object_path_)
    , master(master_)
    , ever_stored(false)
    {
//        cerr << "Creating id " << object_path_ << endl;
    }
    
    boost::shared_ptr<BaseStorage<Object> > clone() const;
            
    void load(boost::shared_ptr<Object> & ptr);
    
    void prefetch(boost::shared_ptr<Object> & ptr) { }
    
    void store(boost::shared_ptr<Object> & ptr);
    
    ~Hdf5Storage();
    
private:
    std::string object_path;
    Hdf5StorageMaster * master;
    bool ever_stored;
    boost::mutex mtx;
};

template<class T, class M>
class locking_ptr
{
public:
    locking_ptr(T * v, M * l)
    : v_(v)
    , l_(l) { }
    
    T * operator->() { return v_; }
    T & operator*() { return *v_; }
    
    ~locking_ptr()
    {
        l_->unlock();
        delete v_;
    }
    
private:
    T * v_;
    M * l_;
        
};

class Hdf5StorageMaster : public BaseStorageMaster
{
public:
    Hdf5StorageMaster(const char * fp)
    : file_path(fp)
    , last_id(0)
    { }
    
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
    
    locking_ptr<alps::hdf5::oarchive, boost::mutex> get_oa_ptr()
    {
        return locking_ptr<alps::hdf5::oarchive, boost::mutex>(new alps::hdf5::oarchive(file_path), &h5mtx);
    }
    
    locking_ptr<alps::hdf5::iarchive, boost::mutex> get_ia_ptr()
    {
        return locking_ptr<alps::hdf5::iarchive, boost::mutex>(new alps::hdf5::iarchive(file_path), &h5mtx);
    }
    
private:
    std::string file_path;
    std::size_t last_id;
    
    boost::mutex h5mtx;
    
    std::set<std::string> prefetch_ids;
    std::set<std::string> load_ids;
    std::set<std::string> write_ids;
    
    boost::condition_variable data_cond;
    boost::mutex data_mut;
};

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
void Hdf5Storage<Object>::load(boost::shared_ptr<Object> & ptr)
{   
//    cerr << "Accessing " << object_path << endl;
    
    static Timer timer("hdf5_load");
    timer.begin();
    
    ptr.reset(new Object());
    locking_ptr<alps::hdf5::iarchive, boost::mutex> ia = master->get_ia_ptr();
    *ia >> alps::make_pvp(object_path, *ptr);
    
    timer.end();
}

template<class Object>
void Hdf5Storage<Object>::store(boost::shared_ptr<Object> & ptr)
{
//    cerr << "Storing " << object_path << endl;
    
    static Timer timer("hdf5_store");
    timer.begin();
    
    locking_ptr<alps::hdf5::oarchive, boost::mutex> oa = master->get_oa_ptr();
    assert(ptr.get() != NULL);
    *oa << alps::make_pvp(object_path, *ptr);
    ptr.reset();
    ever_stored = true;
    
    timer.end();
}

template<class Object>
Hdf5Storage<Object>::~Hdf5Storage()
{
    locking_ptr<alps::hdf5::oarchive, boost::mutex> oa = master->get_oa_ptr();
    if (ever_stored && oa->is_group(object_path))
        oa->delete_group(object_path);
    
//    cerr << "Deleted " << object_path << endl;
}

#endif
