#include "ambient/ambient.h"
#include "ambient/packets/types.h"
#include "ambient/packets/packet.h"
#include "ambient/packets/auxiliary.hpp"
#include "ambient/groups/packet_manager.h"
#include "ambient/groups/group.h"
#include "ambient/core/layout.h"
#include "ambient/auxiliary.hpp"

#include "utils/timings.h"
#include "ambient/core/operation/operation.h"
#include "ambient/core/operation/operation.pp.sa.hpp"

#define AMBIENT_MASTER_RANK 0

using namespace ambient::packets; 
using namespace ambient::groups; 

namespace ambient
{
// global objects accessible anywhere //
    scheduler& layout        = scheduler::instance();
    scheduler& engine        = scheduler::instance();
    multirank& rank          = multirank::instance();
    hash_map&  p_profile_map = hash_map::instance();
    comm_map&  scope_map     = comm_map::instance();
    scope_context& scope     = scope_context::instance();
    access_marker& access    = access_marker::instance();
// global objects accessible anywhere //

    scheduler & scheduler::operator>>(dim2 mem_dim) 
    {
        this->mem_dim  = mem_dim;
        this->default_data_packet_t = new block_packet_t(this->mem_dim*this->item_dim); // to redo in future?
        this->default_data_packet_t->commit();
        if(!world()->get_manager()->subscribed(*this->default_data_packet_t)){
            world()->get_manager()->subscribe(*this->default_data_packet_t);
            world()->get_manager()->add_handler(*this->default_data_packet_t, new core::operation(accept_block, 
                world()->get_manager()->get_pipe(*this->default_data_packet_t, packet_manager::IN)) );
        }
        this->work_dim = NULL;
        this->gpu_dim  = NULL;
        return *this;
    }
    scheduler & scheduler::operator,(dim2 dim) 
    {
        if(this->work_dim == NULL)
            this->work_dim = dim;
        else if(this->gpu_dim == NULL)
            this->gpu_dim = dim;
        return *this;
    }
    scheduler& operator>>(scheduler* instance, dim2 mem_dim) 
    {
        return *instance >> mem_dim;
    }
    scheduler& scheduler::instance()
    {
        static scheduler* singleton = NULL;
        if(!singleton) singleton = new scheduler();
        return *singleton;
    }
    void init()      { engine.init();            }
    void finalize()  { engine.finalize();        }
    void playout()   { engine.playout();         }
    void spin()      { engine.spin();            }
    void spin_loop() { engine.spin_loop();       }
    void world_loop(){ engine.world_loop();      }
    void world_spin(){ engine.world_spin();      }
    int  size()      { return engine.size;       }
    bool occupied()  { return engine.occupied(); }

    scheduler::scheduler(): item_dim(dim2(2,2)), stirring(false){ } // to revert to 128,128
    dim2 scheduler::get_mem_dim() { return this->mem_dim;  }
    dim2 scheduler::get_item_dim(){ return this->item_dim; }
    dim2 scheduler::get_work_dim(){ return this->work_dim; }
    dim2 scheduler::get_gpu_dim() { return this->gpu_dim;  }

    void scheduler::init()
    {
        int threading_level;
        MPI_Init_thread(0, NULL, MPI_THREAD_MULTIPLE, &threading_level);
        MPI_Comm_size(MPI_COMM_WORLD, &this->size);

        this->ambient = new group("ambient", AMBIENT_MASTER_RANK, MPI_COMM_WORLD);
        this->default_data_packet_t = NULL;
        srand48((unsigned int) rank()); // seeding (to redo)
// AUTO TUNING SHOULD START BELOW...

////////////////////////////////////
    }
    void scheduler::finalize()
    {
        MPI_Barrier(this->ambient->mpi_comm);
        MPI_Finalize();
    }
    void scheduler::push(core::operation* logistics, core::operation* computing)
    {
        this->stack.push_back(std::pair<core::operation*,core::operation*>(logistics,computing));
    }
    void scheduler::spin()
    {
        this->world_spin();
        while(!this->router.alt_end_reached()){ // alt is for embedding
            (*this->router.alt_pick())->spin(); // router is slow as it contains hundreds of the same groups
        }
    }
    void scheduler::spin_loop()
    {
        this->world_loop();
        while(!this->router.end_reached()){
            packet_manager* mgr = (*this->router.pick());
            mgr->spin_loop();
        }
    }
    void scheduler::world_spin()
    {
        world()->spin();
    }
    void scheduler::world_loop()
    {
        world()->spin_loop();
    }
    bool scheduler::occupied()
    {
        return this->stirring;
    }
    void scheduler::playout()
    {
        if(this->stack.empty()) return; // easy out
        assert(this->occupied() == false);
        this->stirring = true;
        printf("R%d: playout...\n", ambient::rank());
        MPI_Barrier(this->ambient->mpi_comm);
        one_touch_stack<core::operation*> cleanup_stack;
        std::pair<core::operation*, core::operation*>* pair;
        core::operation* logistics;
        core::operation* computing;
        core::operation* needle_op;
        core::operation* haystack_op;
        bool repeat = true;
        //Timer ct("computing"); Timer dt("distribution"); Timer rt("relation scanning");
        //Timer ambop("ambient pure computing");
        //rt.begin();
        while(!this->stack.end_reached()){
            pair = this->stack.pick();
            pair->first->extract_profiles();
            pair->second->extract_profiles();
        }
// now let's iterate through the stack and mark dependencies
        while(!this->stack.end_reached()){
            needle_op = this->stack.pick()->first;
            do{ haystack_op = this->stack.alt_pick()->first; 
            } while(needle_op != haystack_op); 
            while(!this->stack.alt_end_reached()){
                haystack_op = this->stack.alt_pick()->first;
                for(int i=0; i < needle_op->count; i++)
                for(int j=0; j < haystack_op->count; j++)
                if(needle_op->profiles[i] == haystack_op->profiles[j]){ // pointers comparison
                    if(needle_op->constness[i] && haystack_op->constness[j]) continue;
                    needle_op->add_dependant(haystack_op);
                    goto double_break;
                }
                double_break: continue;
            }
        }
        //rt.end();
// now we all set with dependencies!
        while(repeat)
        {   repeat = false;
            //dt.begin();
            while(!this->stack.end_reached()){
                logistics = this->stack.pick()->first;
                if(logistics->executed) continue;
                if(logistics->dependency_count){ repeat = true; continue; }
                logistics->perform();
                core::apply_changes(logistics->profiles, logistics->count);
                logistics->postprocess();
                if(logistics->get_scope()->involved()){
                    this->router.push_back(logistics->get_scope()->get_manager());
                }
            }
            this->world_loop();
            //dt.end();
            //printf("Computing!\n");
            //ct.begin();
            while(!this->stack.end_reached()){
                pair = this->stack.pick();
                logistics = pair->first;
                computing = pair->second;
                if(logistics->dependency_count || computing->executed) continue;
                cleanup_stack.push_back(logistics);
                if(logistics->get_scope()->involved()){
                    computing->set_scope(logistics->get_scope());
                    scope.set_group(logistics->get_scope());
                    if(logistics->pin == NULL){ // nothing has been pinned
                        logistics->get_scope()->spin_loop();
                        computing->invoke();    // scalapack style
                        logistics->get_scope()->spin_loop();
                    }else{
                        this->context.bind(logistics);
                        //ambop.begin();
                        this->context.discharge(computing);
                        //ambop.end();
                        this->context.finalize();
                    }
                }
                computing->release();
                logistics->release();
            }
            //ct.end();
// cleaning the layout
            while(!cleanup_stack.end_reached()){
                logistics = *cleanup_stack.pick();
                for(int i=0; i < logistics->count; i++)
                    logistics->profiles[i]->clean();
                logistics->resolve_dependencies();
            }
            cleanup_stack.clean(); // :))
            this->router.clean();
        }
        this->stack.clean();
        scope.reset_group();
        this->stirring = false;
    }

    bool is_master(){
        if(scope.get_group() != NULL) return ambient::rank.is_master(scope.get_group());
        return ambient::rank.is_master("ambient");
    }
    group* world(){
        return engine.ambient;
    }

}
