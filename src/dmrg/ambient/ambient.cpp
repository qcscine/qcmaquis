#include "ambient/ambient.h"
#include "ambient/packets/types.h"
#include "ambient/packets/packet.h"
#include "ambient/packets/auxiliary.hpp"
#include "ambient/groups/packet_manager.h"
#include "ambient/groups/group.h"
#include "ambient/groups/auxiliary.hpp"

#define AMBIENT_MASTER_RANK 0

using namespace ambient::packets; 
using namespace ambient::groups; 

namespace ambient
{
// global objects accessible anywhere //
    scheduler& layout = scheduler::instance();
    scheduler& engine = scheduler::instance();
    multirank& rank   = multirank::instance();
    hash_map& void_pt_map = hash_map::instance();
    smp& asmp = smp::instance(); // charge of processes inside kernels

// global objects accessible anywhere //

    scheduler & scheduler::operator>>(dim3 dim_distr) 
    {
        this->dim_distr = dim_distr;
        this->dim_group = NULL;
        this->dim_gpu = NULL;
        return *this;
    }
    scheduler & scheduler::operator,(dim3 dim) 
    {
        if(this->dim_group == NULL){
            this->dim_group = dim;
        }else if(this->dim_gpu == NULL){
            this->dim_gpu = dim;
        }
        return *this;
    }
    scheduler& operator>>(scheduler* instance, dim3 dim_distr) 
    {
        return *instance >> dim_distr;
    }
    scheduler& scheduler::instance()
    {
        static scheduler* singleton = NULL;
        if(!singleton) singleton = new scheduler();
        return *singleton;
    }
    void playout()
    {
        engine.playout(); 
    }
    scheduler::scheduler(): dim_item(dim3(128,128,1))
    {
    }
    dim3 scheduler::group_dim()
    {
        return this->dim_group;
    }
    dim3 scheduler::item_dim()
    {
        return this->dim_item;
    }
    size_t get_bound()
    {
        return 200; // to be redo to something normal
    }
    size_t get_block_bound()
    {
        return 200; // to be redo to something normal
    }
    void scheduler::regression_test()
    {
        printf("Initializing packet system...\n");
        change_t<control_packet_t>(4,3);
        commit_t<control_packet_t>();
        commit_t<data_packet_t>();
        packet* init_packet;
        void* buffer = alloc_t<control_packet_t>();
        if(rank("ambient") == AMBIENT_MASTER_RANK){
            init_packet = pack<control_packet_t>(buffer, 1, "P", rank("ambient"), "DATA", 1);
            send(init_packet, "ambient");
        }else{
            init_packet = recv<control_packet_t>("ambient", buffer);
            printf("%d: Init packet contents: %c %d %c %d %s %d\n", rank("ambient"), init_packet->get<char>(0), 
                                                                    init_packet->get<int>(1), 
                                                                    init_packet->get<char>(2),
                                                                    init_packet->get<int>(3), 
                                                                    (char*)init_packet->get(4), 
                                                                    init_packet->get<int>(5));
        }
        MPI_Barrier(this->comm);

        group* work_grp = new group("work", AMBIENT_MASTER_RANK, this->ambient);
        work_grp->add_every(1);
        work_grp->commit();
        printf("Rank inside work: %d; ambient: %d\n", rank("work"), rank("ambient"));
        int new_ranks[] = { 1, 0 };
        work_grp->reorder(new_ranks);
        work_grp->commit();
        printf("Reordered: Rank inside work: %d; ambient: %d\n", rank("work"), rank("ambient"));
        if(rank("work") == AMBIENT_MASTER_RANK){
            init_packet = pack<control_packet_t>(buffer, 1, "L", rank("work"), "SNCF", 2);
            send(init_packet, "work");
        }else{
            init_packet = recv<control_packet_t>("work", buffer);
            printf("%d: Init packet contents: %c %d %c %d %s %d\n", rank("ambient"), init_packet->get<char>(0), 
                                                                    init_packet->get<int>(1), 
                                                                    init_packet->get<char>(2),
                                                                    init_packet->get<int>(3), 
                                                                    (char*)init_packet->get(4), 
                                                                    init_packet->get<int>(5));
        }
        MPI_Barrier(this->comm);
    }
    void scheduler::init(MPI_Comm comm)
    {
        int threading_level;
        this->comm = comm;
        if(this->comm == NULL){
            MPI_Init_thread(0, NULL, MPI_THREAD_MULTIPLE, &threading_level);
            this->comm = MPI_COMM_WORLD;
        }
        MPI_Comm_size(this->comm, &this->size);

        this->ambient = new group("ambient", AMBIENT_MASTER_RANK, this->comm);

//      regression_test();
        commit_t<control_packet_t>();
        commit_t<data_packet_t>();

/* Distribution kernel example:

void distribution_1(i_dense_matrix* matrix)
{
    group* grp = select(" * from ambient");
    int count = matrix->block_rows*matrix->block_cols / grp->size;
    for(int rank=0; rank < grp->size; rank++){
        for(int j=0; j < count; j++)
        workload[rank] += matrix(rank*count + j); // adding work-group
    }
}

void computation_1(workgroup* block)
{
    block.item[1] += workgroup[10].item[9];
    block += workgroup[7];
}

*/
// AUTO TUNING SHOULD START BELOW...

////////////////////////////////////
    }
    void scheduler::finalize()
    {
        MPI_Barrier(this->comm);
        MPI_Finalize();
    }
    void scheduler::push(core::operation* logistics, core::operation* computing)
    {
        this->logistics_stack.push(logistics);
        this->computing_stack.push(computing);
    }
    void scheduler::playout()
    {
        while(!this->logistics_stack.empty()){
            try{
                this->logistics_stack.front()->perform();
                this->logistics_stack.pop();
            }catch(core::out_of_scope_e e){
                this->logistics_stack.pop();
            }
        }
        while(!this->computing_stack.empty()){
            try{
                this->computing_stack.front()->perform();
                this->computing_stack.pop();
            }catch(core::out_of_scope_e e){
                this->computing_stack.pop();
            }
        }
//        printf("Performing actual communications/computations\n");
    }

    int size(){
        return engine.size;
    }

    bool is_master(){
        return ambient::rank.is_master();
    }
}
