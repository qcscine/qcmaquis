/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *               2011-2012 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#define BOOST_FILESYSTEM_VERSION 3
#include <boost/filesystem.hpp>
#include "stream_storage.h"


/********************
 * StreamStorage
 ********************/
StreamStorage::StreamStorage(std::string const & object_path_,
				             StreamStorageMaster * m)
: object_path(object_path_)
, master_(m)
, status_(Complete)
{ }

StreamStorage::StreamStorage(StreamStorage const & rhs)
: object_path(rhs.master_->get_path())
, master_(rhs.master_)
, status_(StreamStorage::Complete)
{ }

StreamStorage::~StreamStorage()
{
    boost::shared_ptr<StreamRequest> req(new StreamDeleteRequest(this));
    
    master_->add_req(req);
    master_->notify();
}

StreamStorage & StreamStorage::operator=(StreamStorage rhs)
{
    std::swap(*this, rhs);
    return *this;
}

void StreamStorage::wait()
{
    for (std::deque<boost::shared_ptr<boost::mutex> >::iterator it = mutexes.begin();
         it != mutexes.end(); ++it)
        boost::lock_guard<boost::mutex> lock(**it);
    mutexes.clear();
    status_ = Complete;
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

StreamStorageMaster::~StreamStorageMaster()
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

std::string StreamStorageMaster::get_path() const
{
    std::ostringstream oss;
    oss << "/" << last_id++;
    return oss.str();
}

void StreamStorageMaster::print_size() const
{
    if (base_path.size() == 0)
        return;
    
    maquis::cout << "Storage directory size: "; std::cout.flush();
    std::ostringstream oss;
    oss << "du -skh " << base_path << " | awk '{print $1}'";
    system(oss.str().c_str());
}

void StreamStorageMaster::sync()
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

StreamStorage StreamStorageMaster::child()
{
    return StreamStorage(get_path(), this);
}

void StreamStorageMaster::add_req(boost::shared_ptr<StreamRequest> const& req)
{
    boost::lock_guard<boost::mutex> lock(deque_mutex);
    requests.push_back(req);
}

void StreamStorageMaster::notify()
{
    worker_cond.notify_one();
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

