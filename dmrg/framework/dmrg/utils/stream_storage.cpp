#define BOOST_FILESYSTEM_VERSION 3
#include <boost/filesystem.hpp>
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
 * StreamStorageMaster
 ********************/
StreamStorageMaster::StreamStorageMaster(std::string fp, bool enable)
: base_path(fp)
, last_id(0)
, active(true)
, worker(this)
, worker_thread(worker)
{
    
    if (base_path.size() != 0 || enable) {
        
        boost::filesystem::path ph(base_path);
        if (base_path.size() == 0)
            ph = boost::filesystem::temp_directory_path();
        ph = boost::filesystem::unique_path(ph.string() + std::string("/storage_temp_%%%%%%%%%%%%"));
        try {
            boost::filesystem::create_directories(ph);
        } catch (...) {
            maquis::cerr << "Error creating temp dir " << ph << ", try with a different 'storagedir' parameter." << std::endl;
            throw;
        }
        base_path = ph.string();
        maquis::cout << "Temporary storage enabled in " << base_path << std::endl;
    }
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

