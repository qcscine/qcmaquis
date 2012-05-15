#include "ambient/ambient.h"
#include "ambient/models/v_model.h"

namespace ambient { namespace models {

    v_model::reduction::reductionq::reductionq(){

    }

    void v_model::reduction::reductionq::push(v_model::layout::entry* e){

    }

    v_model::reduction::reduction(v_model::revision* r)
    : revision(r)
    {
        dim2 dim = this->revision->get_layout().get_grid_dim();
        for(size_t i = 0; i < dim.y; i++){
            entries.push_back(std::vector<reductionq*>());
            for(size_t j = 0; j < dim.x; j++){
                entries[i].push_back(NULL);
            }
        }
    };

    v_model::reduction::~reduction(){
    }

    v_model::layout::entry& v_model::reduction::operator()(size_t i, size_t j){
        if(this->entries[i][j] == NULL) this->entries[i][j] = new reductionq();
        models::v_model::layout::entry* block = controller.alloc_block(*this->revision);
        this->entries[i][j]->push(block);
        return *block;
    }

} }

