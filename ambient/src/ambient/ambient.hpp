#ifndef AMBIENT_INTERFACE
#define AMBIENT_INTERFACE
#include "ambient/ambient.h"
#include "ambient/utils/memory.hpp"
#include "ambient/utils/jstrings.h"
#include "ambient/interface/forwarding.h"
#include "ambient/utils/io.hpp"
#include "ambient/utils/timings.h"

namespace ambient{

    inline void set_num_threads(size_t n); 

    inline size_t get_num_threads();

    inline void playout();

    inline bool verbose();

    template<typename T>
    inline void assign(const T& ref, int i, int j = 0);

    template<typename T>
    inline void pin(const T& ref, int i, int j = 0);

    template<typename T>
    inline std::pair<size_t*,size_t> id(T& ref);

    template<typename T>
    inline dim2 get_dim(T& ref);

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

    template<typename T>
    inline void assign(const T& ref, int i, int j){ // work_dim = mem_dim
        ambient::models::imodel::revision& revision = current(ref);
        ambient::models::imodel::modifier* op = ctxt.get_op();
        if(op->get_pin() == NULL){
            op->add_condition();
            revision.block(i,j)->get_assignments().push_back(op);
        }
        controller.ifetch_block(revision, i, j);
    }

    template<typename T>
    inline void pin(const T& ref, int i, int j){ // work_dim = mem_dim
        ambient::models::imodel::revision& revision = current(ref);
        ambient::models::imodel::modifier* op = ctxt.get_op();
        op->add_condition();
        revision.block(i,j)->get_assignments().push_back(op);
        controller.ifetch_block(revision, i, j);
    }
    template<ambient::models::imodel::revision&(*STATE)(const ambient::models::v_model::object&), typename T>
    inline void pin(const T& ref, int i, int j){ // work_dim = mem_dim
        ambient::models::imodel::revision& revision = STATE(ref);
        ambient::models::imodel::modifier* op = ctxt.get_op();
        op->add_condition();
        revision.block(i,j)->get_assignments().push_back(op);
        controller.ifetch_block(revision, i, j);
    }

    inline void set_num_threads(size_t n){ 
        ambient::controller.set_num_threads(n);
    }

    inline size_t get_num_threads(){
        return ambient::controller.get_num_threads();
    }

    inline void playout(){
        ambient::controller.flush(); 
    }

    inline bool verbose(){
        return (rank() ? false : true); 
    }

    template<typename T>
    inline std::pair<size_t*,size_t> id(T& ref){
        return current(ref).id();
    }

    template<typename T>
    inline dim2 get_dim(T& ref){
        return current(ref).get_layout().get_dim();
    }
    template<ambient::models::imodel::revision&(*STATE)(const ambient::models::v_model::object&), typename T>
    inline dim2 get_dim(T& ref){
        return STATE(ref).get_layout().get_dim();
    }

    template<typename T>
    inline dim2 get_grid_dim(T& ref){
        return current(ref).get_layout().get_grid_dim();
    }
    template<ambient::models::imodel::revision&(*STATE)(const ambient::models::v_model::object&), typename T>
    inline dim2 get_grid_dim(T& ref){
        return STATE(ref).get_layout().get_grid_dim();
    }

    template<typename T>
    inline dim2 get_mem_grid_dim(T& ref){
        return current(ref).get_layout().get_mem_grid_dim();
    }
    template<ambient::models::imodel::revision&(*STATE)(const ambient::models::v_model::object&), typename T>
    inline dim2 get_mem_grid_dim(T& ref){
        return STATE(ref).get_layout().get_mem_grid_dim();
    }

    template<typename T>
    inline dim2 get_mem_dim(T& ref){
        return current(ref).get_layout().get_mem_dim();
    }
    template<ambient::models::imodel::revision&(*STATE)(const ambient::models::v_model::object&), typename T>
    inline dim2 get_mem_dim(T& ref){
        return STATE(ref).get_layout().get_mem_dim();
    }

    template<typename T>
    inline dim2 get_item_dim(T& ref){
        return current(ref).get_layout().get_item_dim();
    }
    template<ambient::models::imodel::revision&(*STATE)(const ambient::models::v_model::object&), typename T>
    inline dim2 get_item_dim(T& ref){
        return STATE(ref).get_layout().get_item_dim();
    }
    // }}}
}
    
#include "ambient/interface/parallel_t.hpp"
#include "ambient/interface/future_t.hpp"
#endif
