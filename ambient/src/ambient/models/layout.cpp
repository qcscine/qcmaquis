#include "ambient/ambient.h"
#include "ambient/models/v_model.h"

#include "ambient/models/operation/operation.h"
#include "ambient/utils/ceil.h"


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
    : master(0), t_size(t_size)
    {
        this->dim = dim;
    }

    void v_model::layout::embed(void* memory, size_t i, size_t j = 0){

    }

    void v_model::layout::add_path(size_t i, size_t j, size_t node){
        this->paths[i][j]->add_node(node);
    }

    v_model::layout::entry* v_model::layout::get(size_t i, size_t j){
        printf("%d: Trying to access %d x %d of %d x %d\n", this->sid, i, j, this->get_mem_grid_dim().y, this->get_mem_grid_dim().x);
        return this->entries[i][j]; // to correct
    }

    void v_model::layout::mesh(){
        dim2 dim = this->get_mem_grid_dim();
        if(this->mesh_dim.x >= dim.x && this->mesh_dim.y >= dim.y) return;
        for(size_t i = 0; i < dim.y; i++){
            if(i >= this->mesh_dim.y) entries.push_back(std::vector<entry*>());
            for(size_t j = 0; j < dim.x; j++){
                if(j >= this->mesh_dim.x || i >= this->mesh_dim.y) 
                    entries[i].push_back(new v_model::layout::entry());
            }
        }
        if(dim.x > this->mesh_dim.x) this->mesh_dim.x = dim.x;
        if(dim.y > this->mesh_dim.y) this->mesh_dim.y = dim.y;
    }

    void v_model::layout::set_revision(imodel::revision* r){
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
               this->t_size;
    }

    v_model::layout& v_model::layout::operator>>(dim2 dim){
        this->mem_dim = mem_dim;
        this->work_dim = NULL;
        this->item_dim = NULL;
        this->mesh();
        return *this;
    }

    v_model::layout& v_model::layout::operator,(dim2 dim){
        if(this->work_dim == NULL)
            this->work_dim = dim;
        else if(this->item_dim == NULL)
            this->item_dim = dim;
        return *this;
    }

    void v_model::layout::add_route(size_t node){
    
    }

    void v_model::layout::set_dim(dim2 dim){
        if(this->get_mem_grid_dim().y < __a_ceil(dim.y / this->mem_dim.y) || 
           this->get_mem_grid_dim().x < __a_ceil(dim.x / this->mem_dim.x)){
            this->marker.mark(this->get_mem_grid_dim().y, this->get_mem_grid_dim().x);
        }
        this->dim = dim;
        this->mesh();
    }

    dim2 v_model::layout::get_dim() const {
        return this->dim;
    }

    dim2 v_model::layout::get_mem_dim() const { 
        return this->mem_dim;
    }

    dim2 v_model::layout::get_work_dim() const {
        return this->work_dim;
    }

    dim2 v_model::layout::get_item_dim() const { 
        return this->item_dim; 
    }

    dim2 v_model::layout::get_grid_dim() const {
        size_t n = __a_ceil(this->dim.x / this->work_dim.x);
        size_t m = __a_ceil(this->dim.y / this->work_dim.y);
        return dim2(n, m);
    }

    dim2 v_model::layout::get_mem_grid_dim() const {
        size_t n = __a_ceil(this->dim.x / this->mem_dim.x);
        size_t m = __a_ceil(this->dim.y / this->mem_dim.y);
        return dim2(n, m);
    }

    // }}}

    // {{{ layout::entry + layout::marker //

    v_model::layout::entry::entry()
    :header(NULL)
    {
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

    v_model::layout::entry::operator double* (){
        return (double*)this->data;
    }

    v_model::layout::entry::operator std::complex<double>* (){
        return (std::complex<double>*)this->data;
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
