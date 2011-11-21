
#include "stream_storage.h"

/********************
 * StreamStorage
 ********************/
StreamStorage::StreamStorage(std::string const & object_path_,
				             StreamStorageMaster * master_)
: object_path(object_path_)
, master(master_)
, status_(Complete)
{ }

StreamStorage::StreamStorage(StreamStorage const & rhs)
: object_path(rhs.master->get_path())
, master(rhs.master)
, status_(StreamStorage::Complete)
{ }

StreamStorage::~StreamStorage()
{
    boost::shared_ptr<StreamRequest> req(new StreamDeleteRequest(this));
    
    {
        boost::lock_guard<boost::mutex> lock(master->deque_mutex);
        master->requests.push_back(req);
    }
    master->notify();
}


/********************
 * StreamWorker
 ********************/
StreamWorker::StreamWorker(StreamStorageMaster * master_)
: master(master_)
{ }

void StreamWorker::operator()()
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

void StreamWorker::execute()
{
    while (true) {
        boost::shared_ptr<StreamRequest> req;
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

