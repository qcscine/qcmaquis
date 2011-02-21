#ifndef TRIVIAL_SERIALIZER_H
#define TRIVIAL_SERIALIZER_H

template<class Object>
class TrivialStorage : public BaseStorage<Object>
{
public:
    boost::shared_ptr<BaseStorage<Object> > clone() const
    {
        return boost::shared_ptr<BaseStorage<Object> >(new TrivialStorage<Object>());
    }
    void load(boost::shared_ptr<Object>&) { }
    void prefetch(boost::shared_ptr<Object>&) { }
    void store(boost::shared_ptr<Object>&) { }
    ~TrivialStorage() { }
};

class TrivialStorageMaster : public BaseStorageMaster
{
public:
    int flag() { return 0; }
    
    template<class Object>
    boost::shared_ptr<BaseStorage<Object> > make()
    {
        return boost::shared_ptr<BaseStorage<Object> >(new TrivialStorage<Object>());
    }
};

template<> struct storage_master_type<0>
{
    typedef TrivialStorageMaster type;
};

#endif
