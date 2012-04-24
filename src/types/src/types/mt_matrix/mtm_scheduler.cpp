
#include "mtm_scheduler.h"

MtmMaster mtm_global_master;


MtmMaster::MtmMaster(int n)
{
    for (int i = 0; i < n; ++i)
        worker_threads.create_thread(MtmWorker(this));
}

MtmMaster::~MtmMaster ()
{
    worker_threads.interrupt_all();
    notify();
    worker_threads.join_all();
}

void MtmMaster::notify()
{
    worker_cond.notify_one();
}

void MtmMaster::push(boost::shared_ptr<MtmRequest> req)
{
    boost::lock_guard<boost::mutex> lock(queue_mutex);
    queue.push_back(req);
    notify();
}



MtmWorker::MtmWorker(MtmMaster * m) : master(m) { }

void MtmWorker::operator()()
{
    while (true) {
        while (true) {
            boost::unique_lock<boost::mutex> lock(master->queue_mutex);
            //            if (!master->active)
            //                return;
            if (!master->queue.empty())
                break;
            master->worker_cond.wait(lock);
            boost::this_thread::interruption_point();
        }
        
        execute();
    }
}

void MtmWorker::execute()
{
    while (true) {
        boost::shared_ptr<MtmRequest> req;
        {
            boost::lock_guard<boost::mutex> lock(master->queue_mutex);
            if (master->queue.empty())
                break;
            //            cerr << "Queue size: " << master->queue.size() << std::endl;
            req = master->queue.front();
            master->queue.pop_front();
        }
        (*req)();
    }
}

