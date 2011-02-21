#ifndef TEMPORARY_STORAGE_H
#define TEMPORARY_STORAGE_H

#include <boost/shared_ptr.hpp>
#include <boost/noncopyable.hpp>

template<class Object>
class BaseStorage : public boost::noncopyable
{
public:
    virtual boost::shared_ptr<BaseStorage<Object> > clone() const = 0;
    virtual void load(boost::shared_ptr<Object>&) = 0;
    virtual void prefetch(boost::shared_ptr<Object>&) = 0;
    virtual void store(boost::shared_ptr<Object>&) = 0;
    virtual ~BaseStorage() { }
};

class BaseStorageMaster
{
public:
    virtual int flag() = 0;
};

template<int i>
class InvalidStorageMaster
{
public:
    int flag() { return -1; }
    
    template<class Object>
    boost::shared_ptr<BaseStorage<Object> > make()
    {
        std::ostringstream oss;
        oss << "Invalid storage master flag " << i << std::endl;
        throw std::runtime_error(oss.str());
        return boost::shared_ptr<BaseStorage<Object> >();
    }
};

template<int i>
struct storage_master_type
{
    typedef InvalidStorageMaster<i> type;
};

#include "trivial_serializer.h"
#include "hdf5_storage.h"

template<class Object>
boost::shared_ptr<BaseStorage<Object> > storage_factory(BaseStorageMaster & bsm)
{
    switch(bsm.flag())
    {
#define CASE(i) case i: return reinterpret_cast<typename storage_master_type<i>::type*>(&bsm)->make<Object>(); break;
            CASE(0);
            CASE(1);
//            CASE(2);
        default:
            throw std::runtime_error("Undefined storage manager flag.");
            return boost::shared_ptr<BaseStorage<Object> >();
    }
}

template<class Object>
class temporary_storage
{
public:
    temporary_storage(BaseStorage<Object> const & s)
    : object_(new Object())
    , storage_(s.clone())
    , in_memory(true) { }
    
    temporary_storage(BaseStorage<Object> const & s,
                      Object const & o)
    : object_(new Object(o))
    , storage_(s.clone())
    , in_memory(true) { }
    
    temporary_storage(temporary_storage const & rhs)
    : object_(new Object(rhs()))
    , storage_(rhs.storage_->clone())
    , in_memory(true)
    { }
    
    temporary_storage operator=(temporary_storage rhs)
    {
        std::swap(*this, rhs);
        return *this;
    }
    
    temporary_storage operator=(Object const & o)
    {
//        cerr << "Assignment" << endl;
        object_.reset(new Object(o));
        in_memory = true;
        return *this;
    }
    
    Object & operator()()
    {
        if (!in_memory)
            load();
        return *object_;
    }
    Object const & operator()() const
    {
        if (!in_memory)
            load();
        return *object_;
    }
    
    // postcondition: object_ is loaded, valid and unlocked, and all prefetch requests are terminated
    void load() const
    {
        storage_->load(object_);
        in_memory = true;
    }
    
    // postcondition: none
    void prefetch() const
    {
        storage_->prefetch(object_);
    }

    // postcondition: none
    void store() const
    {
        storage_->store(object_);
        in_memory = false;
    }
        
    void swap_with(temporary_storage & rhs)
    {
        swap(this->object_, rhs.object_);
        swap(this->storage_, rhs.storage_);
        swap(this->in_memory, rhs.in_memory);
    }
    
private:
    // const-ness is sort of a joke
    mutable boost::shared_ptr<Object> object_;
    mutable boost::shared_ptr<BaseStorage<Object> > storage_;
    mutable bool in_memory;
};

template<class Object>
void swap(temporary_storage<Object> & a,
          temporary_storage<Object> & b)
{
    a.swap_with(b);
}

#endif
