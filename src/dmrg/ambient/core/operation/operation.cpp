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
    void operation::set_scope(groups::group* scope)
    {
        this->scope = scope;
    }
    void operation::preprocess()
    {
        for(size_t i=0; i < this->count; i++)
            this->profiles[i]->preprocess(ambient::scope.get_group());
        this->set_scope(ambient::scope.get_group());
    }
    groups::group* operation::get_scope()
    {
        return this->scope;
    }
} }
