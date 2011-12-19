#include "ambient/ambient.h"
#include "ambient/core/p_object.h"
#include "ambient/packets/auxiliary.hpp"
#include "ambient/groups/packet_manager.h"

#include "ambient/core/operation/operation.h"
#include "ambient/core/operation/operation.pp.sa.hpp"
#include "ambient/core/auxiliary.h"
    
namespace ambient {

    void accept_block(groups::packet_manager::typed_q& in_q){
        ambient::packets::packet* pack = in_q.get_target_packet();
        try{
            p_object* profile = ambient::model.get_object((unsigned int*)pack->get(A_BLOCK_P_GID_FIELD), 1, pack->get<int>(A_BLOCK_P_ID_FIELD));
            if(pack->get<char>(A_BLOCK_P_STATE_FIELD) == 'P'){
                if(profile->associated_proxy == NULL) throw core::race_condition_e();
                memblock* grp = profile->block(pack->get<int>(A_BLOCK_P_I_FIELD), pack->get<int>(A_BLOCK_P_J_FIELD));
                if(grp->header == NULL) throw core::race_condition_e();
                profile->associated_proxy->reduce( grp, (void*)((size_t)pack->data + profile->get_bound())); // can be done differently
            }else{
                profile->block(pack->get<int>(A_BLOCK_P_I_FIELD), pack->get<int>(A_BLOCK_P_J_FIELD))->set_memory(pack->data);
                pack->lifetime = -1;  // keep it
            }
        }catch(core::race_condition_e){
            pack->lifetime = 2;       // 2 hops because of rethrowing
            in_q.manager->emit(pack); // re-throwing the packet for future handling
        }
    }

    p_object::p_object()
    : dim(NULL), reserved_x(0), reserved_y(0), group_id(0), id(0), init(NULL), default_block(NULL), 
      reference(this), state(ABSTRACT), scope(NULL), xscope(NULL), consted(false), timestamp(0), associated_proxy(NULL), layout(NULL) {
        this->packet_type = ambient::layout.default_data_packet_t;
        this->mem_dim  = engine.get_mem_dim();
        this->item_dim = engine.get_item_dim();
        this->work_dim = engine.get_work_dim();
        this->gpu_dim  = engine.get_gpu_dim();
    };

    p_object::~p_object(){
        if(this->state == SERIAL) return;
        for(int i=0; i < this->skeleton.size(); i++)
            for(int j=0; j < this->skeleton[i].size(); j++)
                delete this->skeleton[i][j];
        delete this->init;
        delete this->layout;
    }

    std::pair<unsigned int*,size_t> p_object::get_id(){
        if(this->id == 0) 
            return std::pair<unsigned int*,size_t>((unsigned int*)&this->id, this->id); // returning valid pointer
        else 
            return std::pair<unsigned int*,size_t>(this->group_id, this->id);
    }

    void p_object::set_scope(groups::group* scope){
        this->xscope = this->scope;
        this->scope = scope;
        this->layout->set_master(scope->get_master_g());
    }
    groups::group* p_object::get_scope(){ return this->scope; }
    groups::group* p_object::get_xscope(){ return this->xscope; } 
    bool p_object::xinvolved(){ return ((this->get_xscope() != NULL && this->get_xscope()->involved()) || this->involved()); } // modified
    bool p_object::involved(){ return this->get_scope()->involved(); }

    p_object & p_object::operator>>(dim2 mem_dim) 
    {
        this->mem_dim = mem_dim;
        this->xpacket_type = this->packet_type;
        this->packet_type = new block_packet_t<typename traits::type>(this->mem_dim*this->item_dim);
        this->packet_type->commit();
        this->reblock();
        this->work_dim = NULL;
        this->gpu_dim = NULL;
        return *this;
    }
    p_object & p_object::operator,(dim2 dim) 
    {
        if(this->work_dim == NULL)
            this->work_dim = dim;
        else if(this->gpu_dim == NULL)
            this->gpu_dim = dim;
        return *this;
    }

    p_object & operator>>(p_object* instance, dim2 mem_dim) {
        return *instance >> mem_dim;
    }

    bool p_object::is_associated_proxy(){
        return (this->state == PROXY);
    }

    void p_object::associate_proxy(void(*R)(memblock*,void*)){
        this->associated_proxy = new p_object();
        this->associated_proxy->reduce    = R;
        this->associated_proxy->reference = this;
        this->associated_proxy->dim       = this->dim;
        this->associated_proxy->mem_dim   = this->mem_dim;
        this->associated_proxy->item_dim  = this->item_dim;
        this->associated_proxy->t_size    = this->t_size; 
        this->associated_proxy->state     = PROXY;
        this->associated_proxy->reblock();
    }

    void p_object::reblock(){
        if(this->state == SERIAL) return;
        int y_size = __a_ceil(this->dim.y / this->get_mem_t_dim().y);
        int x_size = __a_ceil(this->dim.x / this->get_mem_t_dim().x);
        if(this->reserved_x >= x_size && this->reserved_y >= y_size) return;
        for(int i = 0; i < y_size; i++){
            if(i >= this->reserved_y) skeleton.push_back(std::vector<memblock*>());
            for(int j = 0; j < x_size; j++){
                if(j >= this->reserved_x || i >= this->reserved_y) 
                    skeleton[i].push_back(new memblock(this->reference, i, j));
            }
        }
        if(x_size > this->reserved_x) this->reserved_x = x_size;
        if(y_size > this->reserved_y) this->reserved_y = y_size;
    }

    size_t p_object::get_block_lda(){
        return this->get_mem_t_dim().y*this->t_size;
    }


    bool matrix_order_predicate(const core::layout_table::entry& e1, const core::layout_table::entry& e2)
    {
        if(e1.j != e2.j) return e1.j < e2.j;
        return e1.i < e2.i;
    }

    void p_object::solidify(std::vector<core::layout_table::entry> entries)
    {
        int jumper = 0;
// let's find the solid_lda
        std::sort(entries.begin(), entries.end(), matrix_order_predicate);
        this->solid_lda = 0;
        for(int k=0; k < entries.size(); k++) if(entries[k].j == entries[0].j) this->solid_lda++;
        this->data = malloc(__a_ceil(std::max(this->layout->segment_count, this->layout->request_count)/this->solid_lda) *
                            this->get_mem_t_dim().x*this->get_mem_t_dim().y*this->t_size*this->solid_lda);
        assert(this->data != NULL);
        char* memory = (char*)this->data;
// let's copy from separate memory blocks to the general one
        for(int k=0; k < entries.size(); k++){
            for(int j=0; j < this->get_mem_t_dim().x; j++){
                memcpy(memory+j*this->solid_lda*this->get_mem_t_dim().y*this->t_size, 
                       &((char*)this->block(entries[k].i, entries[k].j)->data)[j*this->get_block_lda()], 
                       this->get_mem_t_dim().y*this->t_size);
            }
            memory += this->get_mem_t_dim().y*this->t_size;
            if(++jumper == this->solid_lda){ jumper = 0; memory += this->get_mem_t_dim().y*this->solid_lda*this->t_size*(this->get_mem_t_dim().x-1); }
        }
    }


    void p_object::disperse(std::vector<core::layout_table::entry> entries)
    { 
        int jumper = 0;
        std::sort(entries.begin(), entries.end(), matrix_order_predicate);
        char* memory = (char*)this->data;
        for(int k=0; k < entries.size(); k++){
            for(int j=0; j < this->get_mem_t_dim().x; j++){
                 memcpy(&((char*)this->block(entries[k].i, entries[k].j)->data)[j*this->get_block_lda()], 
                        memory+j*this->solid_lda*this->get_mem_t_dim().y*this->t_size, 
                        this->get_mem_t_dim().y*this->t_size);
            }
            memory += this->get_mem_t_dim().y*this->t_size;
            if(++jumper == this->solid_lda){ jumper = 0; memory += this->get_mem_t_dim().y*this->solid_lda*this->t_size*(this->get_mem_t_dim().x-1); }
        }
        free(this->data);
    }

    void p_object::set_default_block(int i, int j)
    {
        if(i == -1) this->default_block = NULL;
        else this->default_block = this->block(i, j);
    }

    dim2 p_object::get_block_id()
    {
        return dim2(this->default_block->j, this->default_block->i);
    }

    void p_object::constant(){ this->consted = true; }
    void p_object::inconstant(){ this->consted = false; }

    void p_object::preprocess(){
        if(this->state == SERIAL) return;
        groups::group* scope = ambient::scope.get_group();
        ambient::model.update_object(scope, this);
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

    void p_object::postprocess(){
        if(this->state == SERIAL) return;
        if(this->layout != NULL)
        if(this->layout->init_marker.active){
            this->layout->init_marker.clear();
        }
    }
    void p_object::postprocess(int i, int j){
        assert(this->block(i,j)->header == NULL); // if fails - try to reuse this memory
        this->block(i,j)->set_memory(alloc_t(*this->packet_type));
        this->set_default_block(i, j);
        this->init->invoke();
    }

    void p_object::finalize(){
        if(this->state == SERIAL) return;
        if(this->associated_proxy != NULL){
            ambient::scope.spin_loop();
            for(int i=0; i < this->layout->segment_count; i++){ // watch out of constness
                this->get_scope()->get_manager()->emit(pack<layout_packet_t>(alloc_t<layout_packet_t>(), // was scope-manager
                                                                             NULL, "BCAST", "REQUEST FOR REDUCTION DATA",
                                                                             *this->group_id, this->id, "PROXY",
                                                                             ambient::rank(this->get_scope()),
                                                                             this->layout->segment[i].i, 
                                                                             this->layout->segment[i].j));
            }
            ambient::scope.spin_loop();
        }
    }

    void p_object::clean(){
        if(this->state == SERIAL) return;
        this->layout->clean();
        delete this->associated_proxy;
        this->associated_proxy = NULL;
    }

    memblock& p_object::operator()(int i, int j){
        if(this->is_associated_proxy()){ // on-touch init for proxy
            if(!this->block(i,j)->available()){
                this->block(i,j)->set_memory(alloc_t(*this->packet_type));
                memset(this->block(i,j)->data,0,this->get_mem_t_dim().x*this->get_mem_t_dim().y*this->t_size);          
            }
        }else if(!this->block(i,j)->available()){
            groups::packet_manager* manager = world()->get_manager(); //this->consted ? world()->get_manager() : this->get_scope()->get_manager();
            manager->emit(pack<layout_packet_t>(alloc_t<layout_packet_t>(), 
                                                this->layout->get_master(), "P2P", 
                                               "INFORM OWNER ABOUT REQUEST",
                                               *this->group_id, this->id, "GENERIC",
                                                ambient::rank(),      // forward target
                                                i, j));
            while(!this->block(i,j)->available()) ambient::spin();  // spin-lock
        }
        return *(this->block(i,j));
    }

    memblock* p_object::block(int i, int j) const {
        int x_size = __a_ceil(this->dim.x / this->get_mem_t_dim().x);
        int y_size = __a_ceil(this->dim.y / this->get_mem_t_dim().y);
        
        if(i >= y_size || j >= x_size){ printf("Warning: accessing block that is out of range (%d %d)\n", i, j); assert(false); }
        return this->skeleton[i][j];
    }

    size_t p_object::get_bound() const {
        assert(this->packet_type != NULL);
        return (size_t)this->packet_type->displacements[A_BLOCK_P_DATA_FIELD];
    }

    void* p_object::get_data(){
        if(this->default_block == NULL) return NULL; // we asked to convert non structuring arg
        return this->default_block->data;            // >_< need to write proper get for blcck's items
    }
    void p_object::set_init(core::operation* op){
        this->init = op;
    }
    core::operation* p_object::get_init() const {
        return this->init;
    }
    dim2 p_object::get_dim() const {
        return this->dim;
    }
    void p_object::set_dim(dim2 dim){
        if(this->layout != NULL){
            if(this->get_grid_dim().y < __a_ceil(dim.y / this->get_mem_t_dim().y) || 
               this->get_grid_dim().x < __a_ceil(dim.x / this->get_mem_t_dim().x)){
                this->layout->init_marker.mark(this->get_grid_dim().y, this->get_grid_dim().x);
            }
        }
        this->dim = dim;
        this->reblock();
        if(this->layout != NULL)
            this->layout->remap();
    }
    dim2 p_object::get_work_dim() const {
        return this->work_dim;
    }
    void p_object::set_work_dim(dim2 dim){
        this->work_dim = dim;
    }
    dim2 p_object::get_gpu_dim() const {
        return this->gpu_dim;
    }
    void p_object::set_gpu_dim(dim2 dim){
        this->gpu_dim = dim;
    }

    dim2 p_object::get_grid_dim() const {
        int x_size = __a_ceil(this->dim.x / this->get_mem_t_dim().x);
        int y_size = __a_ceil(this->dim.y / this->get_mem_t_dim().y);
        return dim2(x_size, y_size);
    }

    dim2 p_object::get_mem_t_dim() const {
        return this->get_mem_dim() *= this->get_item_dim();
    }
    dim2 p_object::get_mem_dim() const {
        return this->mem_dim;
    }
    void p_object::set_mem_dim(dim2 dim){
        this->mem_dim = dim;
    }

    dim2 p_object::get_item_dim() const {
        return this->item_dim;
    }
    void p_object::set_item_dim(dim2 dim){
        this->item_dim = dim;
    }
}
