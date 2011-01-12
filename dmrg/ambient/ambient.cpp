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

        packet* test_packet;
        if(rank == 0){
//            printf("Size of the type is %d, field 4 is %d\n", (int)sizeof_t<control_packet_t>(), (int)sizeof_t<control_packet_t>(4));
            change_t<control_packet_t>(4, 3);
            commit_t<control_packet_t>();
//            printf("Size of the type is %d, field 4 is %d\n", (int)sizeof_t<control_packet_t>(), (int)sizeof_t<control_packet_t>(4));

            void* buffer = (void*)malloc(sizeof_t<control_packet_t>());
            test_packet = pack<control_packet_t>(buffer, 1, "P", this->rank, "DATA", 1);
            test_packet->send();

        }else{
            change_t<control_packet_t>(4, 3);
            commit_t<control_packet_t>();
            void* buffer = (void*)malloc(sizeof_t<control_packet_t>());

            MPI_Recv(buffer, 1, get_mpi_t<control_packet_t>(), 0, get_t<control_packet_t>()->t_code, this->comm, MPI_STATUS_IGNORE);
            test_packet = unpack<control_packet_t>(buffer);


            printf("Test packet contents %c %d %c %d %s %d\n", *(char*)test_packet->get(0), *(int*)test_packet->get(1),*(char*)test_packet->get(2),*(int*)test_packet->get(3), (char*)test_packet->get(4), *(int*)test_packet->get(5));
        }

        printf("R%d has been initialized\n", this->rank);


// AUTO TUNING SHOULD START BELOW...

////////////////////////////////////
    }

    void scheduler::finalize()
    {
        MPI_Barrier(this->comm);
        MPI_Finalize();
    }
}
