#include "ambient/ambient.h"
#include "ambient/core/operation/operation.h"

namespace ambient{ namespace core{

    void operation::init()
    {
        this->scope = NULL;
        this->pin = NULL;
        this->is_extracted = false;
    }
    void operation::perform()
    {
        ambient::scope.set_op(this);
        ambient::scope.set_group((groups::group*)NULL);
        (this->*prototype)();
        if(!ambient::scope.involved()) return;
        for(size_t i=0; i < this->arg_count; i++)
            this->profiles[i]->postprocess();
    }
    void operation::performx()
    {
        (this->*prototype)();
    }
    void operation::set_ids()
    {
        for(size_t i=0; i < this->arg_count; i++){
            if(this->profiles[i]->id == 0)
                this->profiles[i]->set_id(ambient::scope.get_group()->id);
        }
    }
    void operation::set_scope(groups::group* scope)
    {
        this->scope = scope;
    }
    groups::group* operation::get_scope()
    {
        return this->scope;
    }
} }
