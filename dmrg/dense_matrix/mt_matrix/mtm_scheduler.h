
#ifndef MTM_SCHEDULER_H
#define MTM_SCHEDULER_H

#include <boost/thread.hpp>
#include <boost/shared_ptr.hpp>

class MtmRequest
{
public:
    virtual void operator()() = 0;
    virtual ~MtmRequest() { }
};


class MtmMaster;

class MtmWorker
{
public:
    MtmWorker(MtmMaster * m = NULL);
    
    inline void operator()();
    inline void execute();
    
private:
    MtmMaster * master;
};


class MtmMaster
{
public:
    MtmMaster(int n = 1);
    
    ~MtmMaster();
    
    void notify();
    
    void push(boost::shared_ptr<MtmRequest> req);
    
private:
    std::list<boost::shared_ptr<MtmRequest> > queue;
    boost::mutex queue_mutex;
    
    boost::condition_variable worker_cond;
    
//    template<typename T> friend class mt_matrix;
    friend class MtmWorker;
    
    boost::thread_group worker_threads;
};


extern MtmMaster mtm_global_master;


#endif
