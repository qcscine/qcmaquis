#include "ambient/ambient.h"
#include "ambient/models/v_model.h"
#include "ambient/utils/ceil.h"

#if __APPLE__ && __MACH__
#include <sched.h>	// for sched_yield()
#define pthread_yield() sched_yield()
#endif

//#define LAYOUT_ACCESS_CHECK

namespace ambient { namespace models {

    // {{{ layout model

    v_model::layout::~layout(){
        for(int i = 0; i < this->entries.size(); i++){
            for(int j = 0; j < this->entries[i].size(); j++){
                delete this->entries[i][j];
            }
        }
    }

    v_model::layout::layout(dim2 dim, size_t t_size)
    : master(0), t_size(t_size), mesh_dim(0,0)
    {
        this->dim = dim;
    }

    void v_model::layout::mark(size_t i, size_t j){
        this->marker.mark(i, j);
    }

    bool v_model::layout::marked(size_t i, size_t j){
        return this->marker.marked(i, j);
    }

    void v_model::layout::embed(void* memory, size_t i, size_t j, size_t bound){
        this->get(i,j)->set_memory(memory, bound);
    }

    v_model::layout::entry* v_model::layout::get(size_t i, size_t j){
#ifdef LAYOUT_ACCESS_CHECK
        if(i >= this->get_grid_dim().y || j >= this->get_grid_dim().x)
        printf("%ld: Trying to access %ld x %ld of %ld x %ld\n", this->sid, i, j, this->get_grid_dim().y, this->get_grid_dim().x);
#endif
        return this->entries[i][j];
    }

    void v_model::layout::mesh(){
        this->grid_dim = dim2(__a_ceil(this->dim.x / this->mem_dim.x), 
                              __a_ceil(this->dim.y / this->mem_dim.y));
        dim2 dim = this->grid_dim;
        if(this->mesh_dim.x >= dim.x && this->mesh_dim.y >= dim.y) return;

        for(size_t i = this->mesh_dim.y; i < dim.y; i++)
            entries.push_back(std::vector<entry*>());
        
        for(size_t i = 0; i < this->mesh_dim.y; i++)
            for(size_t j = this->mesh_dim.x; j < dim.x; j++)
                entries[i].push_back(new v_model::layout::entry());

        for(size_t i = this->mesh_dim.y; i < dim.y; i++)
            for(size_t j = 0; j < dim.x; j++)
                entries[i].push_back(new v_model::layout::entry());

        if(dim.x > this->mesh_dim.x) this->mesh_dim.x = dim.x;
        if(dim.y > this->mesh_dim.y) this->mesh_dim.y = dim.y;
    }

    void v_model::layout::set_revision(v_model::revision* r){
        this->revision = static_cast<v_model::revision*>(r);
    }

    std::pair<size_t*,size_t> v_model::layout::id(){
        return std::pair<size_t*,size_t>(this->gid, this->sid);
    }

    size_t v_model::layout::get_master(){
        return this->master;
    }

    size_t v_model::layout::get_mem_size() const {
        return this->get_mem_dim().x *
               this->get_mem_dim().y *
               this->t_size; // returning in bytes
    }
    size_t v_model::layout::get_mem_lda() const {
        return this->get_mem_dim().y *
               this->t_size; // returning lda in bytes
    } 

    v_model::layout& v_model::layout::operator>>(dim2 dim){
        this->mem_dim  = dim;
        this->item_dim = NULL;
        return *this;
    }

    v_model::layout& v_model::layout::operator,(dim2 dim){
        if(this->item_dim == NULL){
            this->item_dim = dim;
            this->mesh();
        }
        return *this;
    }

    void v_model::layout::set_dim(dim2 dim){
        if(this->get_grid_dim().y < __a_ceil(dim.y / this->mem_dim.y) || 
           this->get_grid_dim().x < __a_ceil(dim.x / this->mem_dim.x)){
            this->marker.mark(this->get_grid_dim().y, this->get_grid_dim().x);
        }
        this->dim = dim;
        this->mesh();
    }

    void v_model::layout::set_dimensions(dim2 mem_dim, dim2 item_dim){
        this->mem_dim = mem_dim;
        this->item_dim = item_dim;
        this->mesh();
    }

    dim2 v_model::layout::get_dim() const {
        return this->dim;
    }

    dim2 v_model::layout::get_mem_dim() const { 
        return this->mem_dim;
    }

    dim2 v_model::layout::get_item_dim() const { 
        return this->item_dim; 
    }

    dim2 v_model::layout::get_grid_dim() const {
        return this->grid_dim;
    }

    // }}}

    // {{{ layout::entry + layout::marker //


    bool v_model::layout::entry::trylock(){
        bool acquired = true;
        pthread_mutex_lock(controller.get_pool_control_mutex());
        if(this->locked == false) this->locked = true;
        else acquired = false;
        pthread_mutex_unlock(controller.get_pool_control_mutex());
        return acquired;
    }

    void v_model::layout::entry::unlock(){
        pthread_mutex_lock(controller.get_pool_control_mutex());
        this->locked = false;
        pthread_mutex_unlock(controller.get_pool_control_mutex());
    }

    v_model::layout::entry::~entry(){
        free(this->header);
    }

    v_model::layout::entry::entry()
    : header(NULL), request(false), locked(false)
    {
    }

    v_model::layout::entry::entry(void* memory, size_t bound)
    : header(memory), request(false)
    {
        this->data = (void*)((size_t)memory + bound);
    }

    void v_model::layout::entry::set_memory(void* memory, size_t bound){
        this->header = memory;
        this->data = (void*)((size_t)memory + bound);
    }

    void* v_model::layout::entry::get_memory(){
        return this->header;
    }

    bool v_model::layout::entry::valid(){
        if(this->header != NULL) return true;
        return false;
    }

    bool v_model::layout::entry::requested(){
        return this->request;
    }

    std::list<models::v_model::modifier*>& v_model::layout::entry::get_assignments(){
        return this->assignments;
    }

    std::list<size_t>& v_model::layout::entry::get_path(){
        return this->path;
    }

    v_model::layout::marker::marker(){
        this->active = true;
        this->imarker = 0;
        this->jmarker = 0;
    }

    bool v_model::layout::marker::marked(size_t i, size_t j){
        if(!this->active) return false;
        if(i >= this->imarker || j >= this->jmarker) return true;
        return false;
    }

    void v_model::layout::marker::mark(size_t i, size_t j){
        this->active = true;
        this->imarker = i;
        this->jmarker = j;
    }

    void v_model::layout::marker::clear(){
        this->active = false;
    }

    // }}}

} }
