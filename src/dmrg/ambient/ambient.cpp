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
    hash_map& p_profile_map = hash_map::instance();
    smp& scope = smp::instance(); // charge of processes inside kernels

// global objects accessible anywhere //

    scheduler & scheduler::operator>>(dim3 distr_dim) 
    {
        this->distr_dim = distr_dim;
        this->group_dim = NULL;
        this->gpu_dim = NULL;
        return *this;
    }
    scheduler & scheduler::operator,(dim3 dim) 
    {
        if(this->group_dim == NULL){
            this->group_dim = dim;
            this->default_data_packet_t = new block_packet_t(this->group_dim*this->item_dim); // to redo in future?
            this->default_data_packet_t->commit();
        }else if(this->gpu_dim == NULL){
            this->gpu_dim = dim;
        }
        return *this;
    }
    scheduler& operator>>(scheduler* instance, dim3 distr_dim) 
    {
        return *instance >> distr_dim;
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
    scheduler::scheduler(): item_dim(dim3(128,128,1))
    {
    }
    dim3 scheduler::get_group_dim()
    {
        return this->group_dim;
    }
    dim3 scheduler::get_item_dim()
    {
        return this->item_dim;
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
        this->default_data_packet_t = NULL;
        commit_t<control_packet_t>();
        commit_t<layout_packet_t>();
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
        this->stack.push_back(std::pair<core::operation*,core::operation*>(logistics,computing));
    }
    void scheduler::playout()
    {
        std::pair<core::operation*, core::operation*>* pair;
        core::operation* logistics;
        core::operation* computing;
        while(!this->stack.end_reached())
            this->stack.pick()->first->perform();

        while(!this->stack.end_reached()){
            logistics = this->stack.pick()->first;
            if(logistics->get_scope()->involved()) ambient::core::apply_change_set(logistics->profiles, logistics->count);
            else ambient::core::perform_forwarding(logistics->profiles, logistics->count);
        }

        while(!this->stack.end_reached()){
            pair = this->stack.pick();
            logistics = pair->first;
            if(logistics->get_scope()->involved()){
                computing = pair->second;
                computing->set_scope(logistics->get_scope());
                if(logistics->pin == NULL){ // nothing has been pinned
                    computing->performx();  // scalapack style
                }else{
// performing computation for every item inside every appointed workgroup
                    int pin_cols = logistics->pin->get_grid_dim().x;
                    int pin_rows = logistics->pin->get_grid_dim().y;
                    for(int j=0; j < pin_cols; j++)
                    for(int i=0; i < pin_rows; i++) // todo - extend onto items
                    if((*logistics->pin->layout)(i,j) != NULL){ // rewrite this todo
                        logistics->pin->set_default_group(i, j);
                        computing->performx();
                    }
                } // logistics->pin->set_default_group(-1); // reset in order to avoid mistakes
            }
        }
// cleaning the layout
        while(!this->stack.end_reached()){
            logistics = this->stack.pick()->first;
            for(int i=0; i < logistics->count; i++)
                logistics->profiles[i]->layout->clean();
        }
        this->stack.clean();
    }

    int size(){
        return engine.size;
    }

    bool is_master(){
        return ambient::rank.is_master();
    }
}
