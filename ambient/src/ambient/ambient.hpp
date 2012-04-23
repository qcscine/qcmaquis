#ifndef AMBIENT_INTERFACE_H
#define AMBIENT_INTERFACE_H
#include "ambient/ambient.h"
#include "ambient/utils/memory.hpp"
#include "ambient/utils/jstrings.hpp"
#include "ambient/interface/forwarding.h"
#include "ambient/utils/memory.hpp"
#include "utils/io.hpp"

namespace ambient{

    void set_num_threads(size_t n);

    size_t get_num_threads();

    void playout();

    template<typename T>
    void assign(const T& ref, int i, int j = 0);

    template<typename T>
    void pin(const T& ref, int i, int j = 0);

    template<typename T>
    inline std::pair<size_t*,size_t> id(T& ref);

    template<typename T>
    inline dim2 get_dim(T& ref);

    template<typename T>
    inline dim2 get_work_dim(T& ref);

    template<typename T>
    inline dim2 get_grid_dim(T& ref);

    template<typename T>
    inline dim2 get_mem_grid_dim(T& ref);

    template<typename T>
    inline dim2 get_mem_dim(T& ref);

    template<typename T>
    inline dim2 get_item_dim(T& ref);

    // {{{ realization of interface functions
    #include "ambient/interface/pp/push.pp.hpp" // all variants of push

    void set_num_threads(size_t n){ 
        controller.set_num_threads(n);
    }

    size_t get_num_threads(){
        return controller.get_num_threads();
    }

    void playout(){ 
        controller.flush(); 
    }

    template<typename T>
    void assign(const T& ref, int i, int j){
        ambient::models::imodel::layout& layout = current(ref).get_layout();
        dim2 work_blocks(layout.get_work_dim().x / layout.get_mem_dim().x,
                         layout.get_work_dim().y / layout.get_mem_dim().y);
        size_t ii = i*work_blocks.y;
        size_t jj = j*work_blocks.x;
    
        for(int i = ii; i < ii+work_blocks.y; i++)
            for(int j = jj; j < jj+work_blocks.x; j++){
                if(ctxt.get_op()->get_pin() == NULL){
                    ctxt.get_op()->add_condition();
                    current(ref).block(i,j)->get_assignments().push_back(ctxt.get_op());
                }
                controller.ifetch_block(current(ref), i, j);
            }
    }

    template<typename T>
    void pin(const T& ref, int i, int j){
        ambient::models::imodel::layout& layout = current(ref).get_layout();
        dim2 work_blocks(layout.get_work_dim().x / layout.get_mem_dim().x,
                         layout.get_work_dim().y / layout.get_mem_dim().y);
        size_t ii = i*work_blocks.y;
        size_t jj = j*work_blocks.x;
    
        for(int i = ii; i < ii+work_blocks.y; i++)
            for(int j = jj; j < jj+work_blocks.x; j++){
                ctxt.get_op()->add_condition();
                current(ref).block(i,j)->get_assignments().push_back(ctxt.get_op());
                controller.ifetch_block(current(ref), i, j);
            }
    }

    template<typename T>
    inline std::pair<size_t*,size_t> id(T& ref){
        return current(ref).id();
    }

    template<typename T>
    inline dim2 get_dim(T& ref){
        return current(ref).get_layout().get_dim();
    }

    template<typename T>
    inline dim2 get_grid_dim(T& ref){
        return current(ref).get_layout().get_grid_dim();
    }

    template<typename T>
    inline dim2 get_mem_grid_dim(T& ref){
        return current(ref).get_layout().get_mem_grid_dim();
    }

    template<typename T>
    inline dim2 get_mem_dim(T& ref){
        return current(ref).get_layout().get_mem_dim();
    }

    template<typename T>
    inline dim2 get_item_dim(T& ref){
        return current(ref).get_layout().get_item_dim();
    }

    template<typename T>
    inline dim2 get_work_dim(T& ref){
        return current(ref).get_layout().get_work_dim();
    }

    // }}}
}
    
#include "ambient/interface/parallel_t.hpp"
#endif
