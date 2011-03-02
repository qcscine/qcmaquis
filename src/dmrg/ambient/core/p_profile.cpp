#include "ambient/ambient.h"
#include "ambient/core/p_profile.h"

namespace ambient {

    p_profile::p_profile()
    : reserved_x(0), reserved_y(0), group_id(0), id(0), init_fp(NULL), group_lda(0), default_group(NULL),
      specific(false), profile(this), valid(true), inited(false), need_init(false), master_relay(std::pair<int,int>(-1,-1)), scope(NULL), xscope(NULL) {
        this->packet_type = ambient::layout.default_data_packet_t;
    };

    p_profile* p_profile::dereference(){
        if(!this->valid) printf("Error: attempting to use invalid profile (object was deleted)\n");
        while((this->profile = this->profile->profile) != this->profile->profile);
        return this->profile; // todo - deallocate proxy objects
    }

    void p_profile::set_id(std::pair<unsigned int*,size_t> group_id){
        this->layout = new core::layout_table(this);
        this->group_id = group_id.first;
        this->id = p_profile_map.insert(group_id.first, group_id.second, this->layout);
    }

    void p_profile::set_master(int master){
        this->master_relay = std::pair<int,int>(this->master_relay.second, master);
    }

    int p_profile::get_master(){
        return this->master_relay.second;
    }

    int p_profile::get_xmaster(){
        return this->master_relay.first;
    }

    void p_profile::set_scope(groups::group* scope){
        this->xscope = this->scope;
        this->scope = scope;
    }
    
    groups::group* p_profile::get_scope(){
        return this->scope;
    }

    groups::group* p_profile::get_xscope(){
        return this->xscope;
    }

    p_profile & p_profile::operator>>(dim3 distr_dim) 
    {
        this->specific = true;
        this->distr_dim = distr_dim;
        this->group_dim = NULL;
        this->gpu_dim = NULL;
        return *this;
    }
    p_profile & p_profile::operator,(dim3 dim) 
    {
        if(this->group_dim == NULL){
            this->group_dim = dim;
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

    void p_profile::regroup(){
        if(!this->proxy){
            int y_size = this->dim.y / (this->get_group_dim().y*this->get_item_dim().y);
            int x_size = this->dim.x / (this->get_group_dim().x*this->get_item_dim().x);
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
        }
    }

    size_t p_profile::get_group_lda(){
        if(this->group_lda == 0)
            return this->get_group_dim().y*this->get_item_dim().y*this->type_size;
        else
            return this->group_lda;
    }

    void p_profile::solidify(){
        int i,j;
        size_t offset = 0;
        this->framework = malloc(ambient::get_bound()                            + 
                                 this->layout->count                             *
                                 this->get_group_dim()*this->get_item_dim()      *
                                 this->type_size);
        void* memory = this->data = (void*)((size_t)this->framework + ambient::get_bound());

// let's find the solid_lda
        this->solid_lda = 0; 
        for(j=0; j < this->get_grid_dim().x; j++){
            for(i=0; i < this->get_grid_dim().y; i++) if((*this->layout)(i,j) != NULL) break;
            if((*this->layout)(i,j) != NULL){
                for(i=0; i < this->get_grid_dim().y; i++)
                    if((*this->layout)(i,j) != NULL) this->solid_lda++;
                break;
            }
        }

// let's copy from separate memory blocks to the general one
        for(j=0; j < this->get_grid_dim().x; j++){
            memory = (void*)((size_t)memory + offset*this->get_group_dim().y*this->get_item_dim().y*this->get_group_dim().x*this->get_item_dim().x*this->type_size);
            offset = 0;
            for(i=0; i < this->get_grid_dim().y; i++){
                if((*this->layout)(i,j) != NULL){
                    void* solid_data = (void*)((size_t)memory+offset*this->get_group_dim().y*this->get_item_dim().y*this->type_size);
                    for(int k=0; k < this->get_group_dim().x*this->get_item_dim().x; k++){
                        memcpy((void*)((size_t)solid_data+k*this->solid_lda*this->get_group_dim().y*this->get_item_dim().y), (void*)((size_t)this->group(i,j)->data + k*this->get_group_lda()), 
                               this->get_group_dim().y*this->get_item_dim().y*this->type_size);
                    }
                    offset++;
                }
            }
        }
    }


    void p_profile::disperse(){ 
        size_t offset = 0;
        void* memory = this->data = (void*)((size_t)this->framework + ambient::get_bound());
        for(int j=0; j < this->get_grid_dim().x; j++){
            memory = (void*)((size_t)memory + offset*this->get_group_dim()*this->get_item_dim()*this->type_size);
            offset = 0;
            for(int i=0; i < this->get_grid_dim().y; i++){
                if((*this->layout)(i,j) != NULL){
                    void* solid_data = (void*)((size_t)memory+offset*this->get_group_dim().y*this->get_item_dim().y*this->type_size);
                    for(int k=0; k < this->get_group_dim().x*this->get_item_dim().x; k++){
                        memcpy((void*)((size_t)this->group(i,j)->data + k*this->get_group_lda()), (void*)((size_t)solid_data+k*this->solid_lda*this->get_group_dim().y*this->get_item_dim().y),
                               this->get_group_dim().y*this->get_item_dim().y*this->type_size);
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
        return dim3(this->default_group->i, this->default_group->j, this->default_group->k);
    }

    bool p_profile::is_inited(){
        return this->inited;
    }

    void p_profile::preprocess(groups::group* scope){
        if(this->id == 0) this->set_id(scope->id);
        this->set_master(scope->get_master_g());
        this->set_scope(scope);
        if(this->inited){
            this->need_init = false;
        }else if(!this->inited){
            this->need_init = true;
            this->inited = true;
        }
    }

    void p_profile::postprocess(){
        this->need_init = false;
#ifndef DSCALAPACK_COMPATIBLE
        for(int j=0; j < this->get_grid_dim().x; j++){
            for(int i=0; i < this->get_grid_dim().y; i++){
                if((*this->layout)(i,j) != NULL){
                    this->group(i,j)->set_memory(malloc(ambient::get_block_bound() + this->get_group_dim().x*this->get_group_dim().y *
                                                                                     this->get_item_dim().x*this->get_item_dim().y   *
                                                                                     this->type_size));
                    this->init_fp(this->group(i,j));
                }
            }
        }
#else
        size_t offset = 0;
        this->scope = malloc(ambient::get_bound()                            + 
                             this->layout->count                             *
                             this->get_group_dim()*this->get_item_dim()      *
                             this->type_size);
        void* memory = this->data = (void*)((size_t)this->scope + ambient::get_bound());

        for(int j=0; j < this->get_grid_dim().x; j++){
            memory = (void*)((size_t)memory + offset*this->get_group_dim().y*this->get_item_dim().y*this->get_group_dim().x*this->get_item_dim().x*this->type_size);
            offset = 0;
            for(int i=0; i < this->get_grid_dim().y; i++){
                if((*this->layout)(i,j) != NULL){
                    this->group(i,j)->data = (void*)((size_t)memory+offset*this->get_group_dim().y*this->get_item_dim().y*this->type_size);
                    this->init_fp(this->group(i,j));
                    offset++;
                }
            }
        }
#endif
    }

    workgroup& p_profile::operator()(int i, int j, int k) const {
        return *(this->group(i, j, k));
    }

    workgroup* p_profile::group(int i, int j, int k) const {
        if(this->proxy){
            assert(false); //return new workgroup(&profile, i, j, k);
        }else{
            int x_size = this->dim.x / (this->get_group_dim().x*this->get_item_dim().x);
            int y_size = this->dim.y / (this->get_group_dim().y*this->get_item_dim().y);
            int z_size = this->dim.z / (this->get_group_dim().z*this->get_item_dim().z);

            if(i >= y_size || j >= x_size || k >= z_size) printf("Warning: accessing group that is out of range\n");
            return this->skeleton[i][j];
        }
    }

    void p_profile::imitate(p_profile* profile){
        this->specific   =  profile->specific;
        this->set_gpu_dim  (profile->get_gpu_dim  ());
        this->set_distr_dim(profile->get_distr_dim());
        this->set_group_dim(profile->get_group_dim());
        this->set_item_dim (profile->get_item_dim ());
        this->regroup();
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
        int x_size = this->dim.x / (this->get_group_dim().x * this->get_item_dim().x);
        int y_size = this->dim.y / (this->get_group_dim().y * this->get_item_dim().y);
        int z_size = this->dim.z / (this->get_group_dim().z * this->get_item_dim().z);
        return dim3(x_size, y_size, z_size);
    }

    dim3 p_profile::get_group_dim() const {
        if(this->specific) return this->group_dim;
        else return engine.get_group_dim();
    }
    void p_profile::set_group_dim(dim3 dim){
        this->group_dim = dim;
    }

    dim3 p_profile::get_item_dim() const {
        return engine.get_item_dim();
    }
    void p_profile::set_item_dim(dim3 dim){
        this->item_dim = dim;
    }
    void p_profile::invalidate(){
        if(!this->proxy) this->valid = false;
    }
    bool p_profile::is_valid(){
        return this->valid;
    }
}
