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
        (this->*prototype)();
    }
    void operation::invoke()
    {
        (this->*prototype)();
    }
    void operation::preprocess()
    {
        this->set_scope(ambient::scope.get_group());
        for(size_t i=0; i < this->count; i++) this->profiles[i]->preprocess();
    }
    void operation::finalize()
    {
        for(size_t i=0; i < this->count; i++) this->profiles[i]->finalize();
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
