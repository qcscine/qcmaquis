#include "ambient/ambient.h"
#include "ambient/types.h"
#include "ambient/packet.h"
#include "ambient/packet_manager.h"

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
//        printf("Size of the type is %d, field 4 is %d\n", (int)packet_type::get_size<control_packet_type>(), (int)packet_type::get_size<control_packet_type>(4));
        packet_type::change_size<control_packet_type>(4, 3);
//        printf("Size of the type is %d, field 4 is %d\n", (int)packet_type::get_size<control_packet_type>(), (int)packet_type::get_size<control_packet_type>(4));

        packet_type::get<control_packet_type>()->commit();
        test_packet = new packet(packet_type::get<control_packet_type>(), 1, "P", this->rank, "DATA", 1);
        test_packet->send();

    }else{
        packet_type::change_size<control_packet_type>(4, 3);
        packet_type::get<control_packet_type>()->commit();
        void* buffer = (void*)malloc(packet_type::get<control_packet_type>()->type_size);

        MPI_Recv(buffer, 1, packet_type::get<control_packet_type>()->mpi_type, 0, packet_type::get<control_packet_type>()->type_code, this->comm, MPI_STATUS_IGNORE);
        test_packet = new packet(buffer);


        printf("Test packet contents %c %d %c %d %s %d\n", *(char*)test_packet->get(0), *(int*)test_packet->get(1),*(char*)test_packet->get(2),*(int*)test_packet->get(3), (char*)test_packet->get(4), *(int*)test_packet->get(5));
    }

//    MPI_Aint extent; int sizeT;
//    MPI_Type_extent(test_packet->get_mpi_type(), &extent);
//    printf("INT extent is %d\n", (int)extent);
//    MPI_Type_size(test_packet->get_mpi_type(), &sizeT);
//    printf("INT size is %d\n", sizeT);
    printf("R%d has been initialized\n", this->rank);


// AUTO TUNING SHOULD START BELOW...

////////////////////////////////////
}

void scheduler::finalize()
{
    MPI_Barrier(this->comm);
    MPI_Finalize();
}
