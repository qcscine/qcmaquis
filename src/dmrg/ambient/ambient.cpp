#include "ambient/ambient.h"
#include "ambient/packets/types.h"
#include "ambient/packets/packet.h"
#include "ambient/packets/auxiliary.h"
#include "ambient/groups/packet_manager.h"
#include "ambient/groups/group.h"
#include "ambient/auxiliary.h"

#define AMBIENT_MASTER_RANK 0

using namespace ambient::packets; 
using namespace ambient::groups; 

namespace ambient
{
    scheduler& scheduler::instance(){
        static scheduler* singleton = NULL;
        if(!singleton) singleton = new scheduler();
        return *singleton;
    }
    scheduler::scheduler():rank(multirank::instance())
    {
    }

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

        this->nest = new group("nest", AMBIENT_MASTER_RANK, this->comm);
        if(this->rank("nest") == AMBIENT_MASTER_RANK) this->mode = AMBIENT_MASTER;
        else this->mode = GROUP_SLAVE;


// initializing our data-types:

        printf("Initializing packet system...\n");
        change_t<control_packet_t>(4,3);
        commit_t<control_packet_t>();
        commit_t<data_packet_t>();
        packet* init_packet;
        void* buffer = alloc_t<control_packet_t>();
        if(this->rank("nest") == AMBIENT_MASTER_RANK){
            init_packet = pack<control_packet_t>(buffer, 1, "P", this->rank("nest"), "DATA", 1);
            send(init_packet, this->nest);
        }else{
            init_packet = recv<control_packet_t>(this->nest, buffer);
            printf("Init packet contents: %c %d %c %d %s %d\n", init_packet->get<char>(0), 
                                                                init_packet->get<int>(1), 
                                                                init_packet->get<char>(2),
                                                                init_packet->get<int>(3), 
                                                                (char*)init_packet->get(4), 
                                                                init_packet->get<int>(5));
        }
        MPI_Barrier(this->comm);

        printf("initialized\n");


// AUTO TUNING SHOULD START BELOW...

////////////////////////////////////
    }

    void scheduler::finalize()
    {
        MPI_Barrier(this->comm);
        MPI_Finalize();
    }
}
