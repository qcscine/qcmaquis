#include "ambient/ambient.h"
#include "ambient/packets/types.h"
#include "ambient/packets/packet.h"
#include "ambient/packets/packet_manager.h"
#include "ambient/packets/auxiliary.h"

using namespace ambient::packets; 

namespace ambient
{
    scheduler* scheduler::instance(){
        static scheduler* singleton = NULL;
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


// initializing our data-types:

        print<0>("Initializing packet system...\n");
        packet_manager::instance()->set_comm(MPI_COMM_WORLD);
        packet* init_packet;
        change_t<control_packet_t>(4,3);
        commit_t<control_packet_t>();
        commit_t<data_packet_t>();
        void* buffer = alloc_t<control_packet_t>();
        if(rank == 0){
            init_packet = pack<control_packet_t>(buffer, 1, "P", this->rank, "DATA", 1);
            init_packet->send();
        }else{
            init_packet = recv<control_packet_t>(buffer);
            printf("Init packet contents: %c %d %c %d %s %d\n", init_packet->get<char>(0), 
                                                                init_packet->get<int>(1), 
                                                                init_packet->get<char>(2),
                                                                init_packet->get<int>(3), 
                                                                (char*)init_packet->get(4), 
                                                                init_packet->get<int>(5));
        }
        MPI_Barrier(this->comm);
        print("initialized\n");


// AUTO TUNING SHOULD START BELOW...

////////////////////////////////////
    }

    void scheduler::finalize()
    {
        MPI_Barrier(this->comm);
        MPI_Finalize();
    }
}
