#include "ambient/ambient.h"
#include "ambient/core/operation/operation.h"
#include "ambient/auxiliary.hpp"

namespace ambient{ namespace core{

    void operation::add_dependant(operation* dep){
        if(this->dependants == NULL) 
            this->dependants = new one_touch_stack<operation*>();
        this->dependants->push_back(dep);
        dep->dependency_count++;
    }

    void operation::resolve_dependencies(){
        if(this->dependants == NULL) return; 
        while(!this->dependants->end_reached()){
            (*this->dependants->pick())->dependency_count--;
        }
    }

    void operation::init()
    {
        one_touch_stack<operation*> deps;
        this->scope = NULL;
        this->pin = NULL;
        this->dependants = NULL;
        this->dependency_count = 0;
        this->executed = false;
    }
    void operation::perform()
    {
        ambient::scope.set_op(this);
        (this->*prototype)();
        this->executed = true;
    }
    void operation::invoke()
    {
        (this->*prototype)();
        this->executed = true;
    }
    void operation::extract_profiles()
    {
        (this->*extract)();
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
