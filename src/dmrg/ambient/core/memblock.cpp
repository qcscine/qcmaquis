#include "ambient/ambient.h"
#include "ambient/core/memblock.h"

namespace ambient {

    memblock::~memblock(){
        //free(this->header);
    }

    memblock::memblock(p_profile* p, int i, int j)
    : profile(p), i(i), j(j), header(NULL), data(NULL), timestamp(0) {};

    p_profile* memblock::get_profile(){
        return this->profile;
    }

    dim2 memblock::get_mem_dim(){
        return this->get_profile()->get_mem_dim();
    }

    dim2 memblock::get_mem_t_dim(){
        return this->get_profile()->get_mem_t_dim();
    }

    dim2 memblock::get_item_dim(){
        return this->get_profile()->get_item_dim();
    }

    bool memblock::available(){
        return (this->timestamp == this->get_profile()->timestamp);
    }

    void memblock::set_memory(void* memory){
        this->header = memory;
        this->data = (void*)((size_t)memory + this->get_profile()->get_bound());
        this->timestamp = this->get_profile()->timestamp;
    }

    void* memblock::element(int i, int j){
        int lda = this->profile->get_mem_t_dim().y;
        return &((char*)this->data)[(j*lda+i)*this->profile->t_size];
    }

    void* memblock::operator()(int i, int j){
        return this->item(i,j);
    }

    void* memblock::item(int i, int j){
        i = i*this->profile->get_item_dim().y;
        j = j*this->profile->get_item_dim().x;
    
        int x_size = this->profile->get_mem_dim().x;
        int y_size = this->profile->get_mem_dim().y;
        if(i >= y_size || j >= x_size) printf("Warning: accessing block item that is out of range (%d %d)\n", i, j);
        return (void*)((size_t)this->data + j*y_size + i); // the model of accessing actual data can be changed in future - will need to try out!
    }
}
