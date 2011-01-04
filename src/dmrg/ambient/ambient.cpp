#include "ambient/ambient.h"

scheduler* scheduler::singleton = NULL;
scheduler* scheduler::instance(){
    if(!singleton) singleton = new scheduler();
    return singleton;
}
scheduler::scheduler(){};

scheduler & scheduler::operator>>(dim3 dim_distr) 
{
    this->dim_distr = dim_distr;
    this->dim_cpu = NULL;
    this->dim_gpu = NULL;
    return *this;
}
scheduler & scheduler::operator,(dim3 dim) 
{
    if(this->dim_cpu == NULL){
        this->dim_cpu = dim;
    }else if(this->dim_gpu == NULL){
        this->dim_gpu = dim;
    }
    return *this;
}

scheduler& operator>>(scheduler* instance, dim3 dim_distr) {
    return *instance >> dim_distr;
}

void scheduler::initialize(MPI_Comm comm)
{
    int threading_level;
    this->comm = comm;
    if(this->comm == NULL){
        MPI_Init_thread(0, NULL, MPI_THREAD_MULTIPLE, &threading_level);
        this->comm = MPI_COMM_WORLD;
    }
    MPI_Comm_size(this->comm, &this->size);
    MPI_Comm_rank(this->comm, &this->rank);

    if(this->rank == 0) this->mode = AMBIENT_MASTER;
    else this->mode = GROUP_SLAVE;

    printf("R%d has been initialized\n", this->rank);

// AUTO TUNING SHOULD START BELOW...

////////////////////////////////////



}

void scheduler::finalize()
{
    MPI_Barrier(this->comm);
    MPI_Finalize();
    delete singleton;
    singleton = NULL;
}
