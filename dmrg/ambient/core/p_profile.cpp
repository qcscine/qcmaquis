#include "ambient/ambient.h"
#include "ambient/core/p_profile.h"
#include "ambient/packets/auxiliary.hpp"
#include "ambient/groups/packet_manager.h"

#include "ambient/core/operation/operation.h"
#include "ambient/core/operation/operation.pp.sa.hpp"
#include "ambient/core/auxiliary.h"

namespace ambient {
    void accept_block(groups::packet_manager::typed_q& in_q){
        ambient::packets::packet* pack = in_q.get_target_packet();
        try{
            p_profile* profile = p_profile_map.find((unsigned int*)pack->get(A_BLOCK_P_GID_FIELD), 1, pack->get<int>(A_BLOCK_P_ID_FIELD))->profile;
            if(pack->get<char>(A_BLOCK_P_STATE_FIELD) == 'P'){
                if(profile->associated_proxy == NULL) throw core::race_condition_e();
                workgroup* grp = profile->group(pack->get<int>(A_BLOCK_P_I_FIELD), pack->get<int>(A_BLOCK_P_J_FIELD));
                if(grp->header == NULL) throw core::race_condition_e();
                profile->associated_proxy->reduce_fp( grp, (void*)((size_t)pack->data + profile->get_bound())); // can be done differently
            }else{
                profile->group(pack->get<int>(A_BLOCK_P_I_FIELD), pack->get<int>(A_BLOCK_P_J_FIELD))->set_memory(pack->data);
            }
        }catch(core::race_condition_e){
            in_q.manager->emit(pack); // re-throwing the packet for future handling
        }
    }

    p_profile::p_profile()
    : reserved_x(0), reserved_y(0), group_id(0), id(0), init_fp(NULL), group_lda(0), default_group(NULL),
      profile(this), valid(true), state(ABSTRACT), master_relay(std::pair<int,int>(-1,-1)), scope(NULL), xscope(NULL), consted(false), timestamp(0), associated_proxy(NULL) {
        this->packet_type = ambient::layout.default_data_packet_t;
        this->group_dim = engine.get_group_dim();
        this->item_dim  = engine.get_item_dim();
        this->distr_dim = engine.get_distr_dim();
        this->gpu_dim   = engine.get_gpu_dim();
    };

    p_profile* p_profile::dereference(){
        if(!this->valid) printf("Error: attempting to use invalid profile (object was deleted)\n");
        while((this->profile = this->profile->profile) != this->profile->profile);
        return this->profile; // todo - deallocate proxy objects
    }

    void p_profile::operator=(const p_profile& profile){
        // todo I guess
        this->profile      = const_cast<p_profile*>(&profile);
        this->state        = PROXY;                   // to handle properly
        this->group_id     = profile.group_id;
        this->id           = profile.id;
        this->scope        = profile.scope;
        this->layout       = profile.layout;         // pointer
        this->dim          = profile.dim;
        this->t_size       = profile.t_size; 
        this->packet_type  = profile.packet_type;    // pointer
        this->distr_dim    = profile.get_distr_dim();
        this->group_dim    = profile.get_group_dim();
        this->item_dim     = profile.get_item_dim(); 
        this->gpu_dim      = profile.get_gpu_dim();  
        //this->skeleton; 

        //this->xscope       = profile.xscope;       // not needed
        //this->data         = profile.data;         // not needed
        //this->master_relay = profile.master_relay; // not needed ?
        //this->timestamp    = profile.timestamp;    // not needed
        //this->lda          = profile.lda;          // not needed
        //this->solid_lda    = profile.solid_lda;    // not needed
        //this->group_lda    = profile.group_lda;    // not needed
        //this->state        = profile.state;        // not needed ?
        //this->reserved_x   = profile.reserved_x;   // not needed
        //this->reserved_y   = profile.reserved_y;   // not needed
    }

    p_profile* p_profile::associate_proxy(p_profile* proxy, void(*R)(workgroup*,void*)){
        this->associated_proxy = proxy;
        this->associated_proxy->reduce_fp = R;
    }

    void p_profile::set_id(std::pair<unsigned int*,size_t> group_id){
        this->layout = new core::layout_table(this);
        this->group_id = group_id.first;
        this->id = p_profile_map.insert(group_id.first, group_id.second, this->layout);
    }

    void p_profile::set_master(int master){
        this->master_relay = std::pair<int,int>(this->master_relay.second, master);
    }

    void p_profile::set_scope(groups::group* scope){
        this->xscope = this->scope;
        this->scope = scope;
        this->set_master(scope->get_master_g());
    }
    
    int p_profile::get_master(){ return this->master_relay.second; }
    int p_profile::get_xmaster(){ return this->master_relay.first; }
    groups::group* p_profile::get_scope(){ return this->scope; }
    groups::group* p_profile::get_xscope(){ return this->xscope; } 
    bool p_profile::xinvolved(){ return ((this->get_xscope() != NULL && this->get_xscope()->involved()) || this->involved()); } // modified
    bool p_profile::involved(){ return this->get_scope()->involved(); }

    p_profile & p_profile::operator>>(dim3 distr_dim) 
    {
        this->distr_dim = distr_dim;
        this->group_dim = NULL;
        this->gpu_dim = NULL;
        return *this;
    }
    p_profile & p_profile::operator,(dim3 dim) 
    {
        if(this->group_dim == NULL){
            this->group_dim = dim;
            this->xpacket_type = this->packet_type;
            this->packet_type = new block_packet_t(this->group_dim*this->item_dim);
            this->packet_type->commit();
            this->regroup();
        }else if(this->gpu_dim == NULL){
            this->gpu_dim = dim;
        }
        return *this;
    }

    p_profile & operator>>(p_profile* instance, dim3 distr_dim) {
        return *instance >> distr_dim;
    }

    bool p_profile::is_proxy(){
        return (this->state == PROXY);
    }

    void p_profile::regroup(){
        if(!this->is_proxy()){
            int y_size = this->dim.y / (this->get_group_t_dim().y);
            int x_size = this->dim.x / (this->get_group_t_dim().x);
            if(this->reserved_x >= x_size && this->reserved_y >= y_size) return;
            for(int i = 0; i < y_size; i++){
                if(i >= this->reserved_y) skeleton.push_back(std::vector<workgroup*>());
                for(int j = 0; j < x_size; j++){
                    if(j >= this->reserved_x || i >= this->reserved_y) 
                        skeleton[i].push_back(new workgroup(&profile, i, j));
                }
            }
            if(x_size > this->reserved_x) this->reserved_x = x_size;
            if(y_size > this->reserved_y) this->reserved_y = y_size;
        }else{
        }
    }

    size_t p_profile::get_group_lda(){
        return this->get_group_t_dim().y*this->t_size;
    }

    void p_profile::solidify(){
        int i,j;
        size_t offset = 0;
        this->data = malloc(this->layout->segment_count                  *
                            (this->get_group_dim()*this->get_item_dim()) *
                            this->t_size);
        void* memory = this->data;

// let's find the solid_lda
        this->solid_lda = 0; 
        for(j=0; j < this->get_grid_dim().x; j++){
            for(i=0; i < this->get_grid_dim().y; i++){
                if(this->group(i,j)->available()){
                    for(i=0; i < this->get_grid_dim().y; i++)
                        if(this->group(i,j)->available()) this->solid_lda++;
                    break;
                }
            }
            if(this->solid_lda) break;
        }
// let's copy from separate memory blocks to the general one
        for(j=0; j < this->get_grid_dim().x; j++){
            memory = (void*)((size_t)memory + offset*(this->get_group_dim()*this->get_item_dim())*this->t_size);
            offset = 0;
            for(i=0; i < this->get_grid_dim().y; i++){
                if(this->group(i,j)->available()){
                    void* solid_data = (void*)((size_t)memory+offset*this->get_group_t_dim().y*this->t_size);
                    for(int k=0; k < this->get_group_t_dim().x; k++){
                        assert(this->group(i,j)->data != NULL);
                        memcpy((void*)((size_t)solid_data+k*this->solid_lda*this->get_group_t_dim().y), (void*)((size_t)this->group(i,j)->data + k*this->get_group_lda()), 
                               this->get_group_t_dim().y*this->t_size);
                    }
                    offset++;
                }
            }
        }
    }


    void p_profile::disperse(){ 
        size_t offset = 0;
        void* memory = this->data;
        for(int j=0; j < this->get_grid_dim().x; j++){
            memory = (void*)((size_t)memory + offset*(this->get_group_dim()*this->get_item_dim())*this->t_size);
            offset = 0;
            for(int i=0; i < this->get_grid_dim().y; i++){
                if(this->group(i,j)->available()){
                    void* solid_data = (void*)((size_t)memory+offset*this->get_group_t_dim().y*this->t_size);
                    for(int k=0; k < this->get_group_t_dim().x; k++){
                        memcpy((void*)((size_t)this->group(i,j)->data + k*this->get_group_lda()), (void*)((size_t)solid_data+k*this->solid_lda*this->get_group_t_dim().y),
                               this->get_group_t_dim().y*this->t_size);
                    }
                    offset++;
                }
            }
        }
    }

    void p_profile::set_default_group(int i, int j, int k)
    {
        if(i == -1) this->default_group = NULL;
        else this->default_group = this->group(i, j, k);
    }

    dim3 p_profile::get_group_id()
    {
        return dim3(this->default_group->j, this->default_group->i, this->default_group->k);
    }

    void p_profile::touch(){
        
        if(this->state == ABSTRACT){
            this->state = COMPOSING;
        }else if(this->state == COMPOSING){
            this->state = GENERIC;
        }

    }

    void p_profile::constant(){ this->consted = true; }
    void p_profile::inconstant(){ this->consted = false; }

    void p_profile::preprocess(){
        groups::group* scope = ambient::scope.get_group();
        if(this->id == 0) this->set_id(scope->id);
        this->touch();
        if(!this->consted || this->state == COMPOSING){ // bad case - the same argument twice :/
            this->timestamp++;
            this->set_scope(scope);
        }
        if(scope->involved()){
            if(!scope->get_manager()->subscribed(*this->packet_type)){
                scope->get_manager()->subscribe(*this->packet_type);
                scope->get_manager()->add_handler(*this->packet_type, new core::operation(accept_block, 
                    scope->get_manager()->get_pipe(*this->packet_type, groups::packet_manager::IN)) );
            }
            if(!world()->get_manager()->subscribed(*this->packet_type)){
                world()->get_manager()->subscribe(*this->packet_type);
                world()->get_manager()->add_handler(*this->packet_type, new core::operation(accept_block, 
                    world()->get_manager()->get_pipe(*this->packet_type, groups::packet_manager::IN)) );
            }
        }
    }

    void p_profile::postprocess(){
        int i, j;
        for(int k=0; k < this->layout->segment_count; k++){
            i = this->layout->segment[k].i;
            j = this->layout->segment[k].j;
            if(this->group(i,j)->header != NULL) continue; // avoiding redunant allocations
            this->group(i,j)->set_memory(alloc_t(*this->packet_type));
            this->init_fp(this->group(i,j));
        }
    }

    void p_profile::finalize(){
        if(this->associated_proxy != NULL){
            for(int i=0; i < this->layout->segment_count; i++){
                this->get_scope()->get_manager()->emit(pack<layout_packet_t>(alloc_t<layout_packet_t>(), // was scope-manager
                                                                             NULL, "BCAST", "REQUEST FOR REDUCTION DATA",
                                                                             *this->group_id, this->id, "PROXY",
                                                                             ambient::rank(this->get_scope()),
                                                                             this->layout->segment[i].i, 
                                                                             this->layout->segment[i].j, 
                                                                             this->layout->segment[i].k));
            }
        }
    }

    void p_profile::clean(){
        this->layout->clean();
        delete this->associated_proxy;
        this->associated_proxy = NULL;
    }

    workgroup& p_profile::operator()(int i, int j, int k){
        if(this->is_proxy()){ // on-touch init for proxy
            this->group(i,j,k)->set_memory(alloc_t(*this->packet_type));
        }else if(!this->group(i,j,k)->available()){
            //printf("R%d: Requesting the group which is outdated (%d: %d %d)\n", ambient::rank(), this->id, i, j); // let's make the request here!
            groups::packet_manager* manager = world()->get_manager(); //this->consted ? world()->get_manager() : this->get_scope()->get_manager();
            manager->emit(pack<layout_packet_t>(alloc_t<layout_packet_t>(), 
                                                this->get_master(), "P2P", 
                                               "INFORM OWNER ABOUT REQUEST",
                                               *this->group_id, this->id, "GENERIC",
                                                ambient::rank(),      // forward target
                                                i, j, k));
            while(!this->group(i,j,k)->available()) ambient::spin();  // spin-lock
        }
        return *(this->group(i, j, k));
    }

    workgroup* p_profile::group(int i, int j, int k) const {
        int x_size = this->dim.x / this->get_group_t_dim().x;
        int y_size = this->dim.y / this->get_group_t_dim().y;
        int z_size = this->dim.z / this->get_group_t_dim().z;
        
        if(i >= y_size || j >= x_size || k >= z_size) printf("Warning: accessing group that is out of range (%d %d %d)\n", i, j, k);
        return this->skeleton[i][j];
    }

    void p_profile::imitate(p_profile* profile){
        this->set_gpu_dim  (profile->get_gpu_dim  ());
        this->set_distr_dim(profile->get_distr_dim());
        this->set_group_dim(profile->get_group_dim());
        this->set_item_dim (profile->get_item_dim ());
        this->regroup();
    }

    size_t p_profile::get_bound() const {
        assert(this->packet_type != NULL);
        return (size_t)this->packet_type->displacements[A_BLOCK_P_DATA_FIELD];
    }

    void* p_profile::get_data(){
        if(this->default_group == NULL) return NULL; // we asked to convert non structuring arg
        return this->default_group->data;            // >_< need to write proper get for group's items
    }
    dim3 p_profile::get_dim() const {
        return this->dim;
    }
    void p_profile::set_dim(dim3 dim){
        this->dim = dim;
    }
    dim3 p_profile::get_distr_dim() const {
        return this->distr_dim;
    }
    void p_profile::set_distr_dim(dim3 dim){
        this->distr_dim = dim;
    }
    dim3 p_profile::get_gpu_dim() const {
        return this->gpu_dim;
    }
    void p_profile::set_gpu_dim(dim3 dim){
        this->gpu_dim = dim;
    }

    dim3 p_profile::get_grid_dim() const {
        int x_size = this->dim.x / this->get_group_t_dim().x;
        int y_size = this->dim.y / this->get_group_t_dim().y;
        int z_size = this->dim.z / this->get_group_t_dim().z;
        return dim3(x_size, y_size, z_size);
    }

    dim3 p_profile::get_group_t_dim() const {
        return this->get_group_dim() *= this->get_item_dim();
    }
    dim3 p_profile::get_group_dim() const {
        return this->group_dim;
    }
    void p_profile::set_group_dim(dim3 dim){
        this->group_dim = dim;
    }

    dim3 p_profile::get_item_dim() const {
        return this->item_dim;
    }
    void p_profile::set_item_dim(dim3 dim){
        this->item_dim = dim;
    }
    void p_profile::invalidate(){
        if(!this->is_proxy()) this->valid = false;
    }
    bool p_profile::is_valid(){
        return this->valid;
    }
}
