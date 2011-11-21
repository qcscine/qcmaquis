#include "ambient/ambient.h"
#include "ambient/core/operation/operation.h"
#include "ambient/auxiliary.hpp"

namespace ambient{ namespace core{

    operation::~operation(){
        free(this->profiles);
        free(this->constness);
        delete[] this->arguments;
        delete this->dependants;
    }

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
        this->constness = NULL;
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
    void operation::extract_profiles()
    {
        (this->*extract)();
        this->constness = (bool*)malloc(sizeof(bool)*count);
        this->extract_constness();
    }
    void operation::extract_constness(){
        for(int i=0; i < this->count; i++){
            this->constness[i] = this->profiles[i]->consted;
        } // playout can be further enhanced for ABSTRACT state
    }
    void operation::preprocess()
    {
        this->set_scope(ambient::scope.get_group());
        for(size_t i=0; i < this->count; i++){
            this->profiles[i]->consted = this->constness[i] && this->profiles[i]->state != ABSTRACT; // regaining const info
            this->profiles[i]->preprocess();
        }
    }
    void operation::postprocess()
    {
        for(size_t i=0; i < this->count; i++)
            this->profiles[i]->postprocess();
    }
    void operation::finalize()
    {
        for(size_t i=0; i < this->count; i++) this->profiles[i]->finalize();
    }
    void operation::release() // called only for one op in a pair
    {
        for(int i=0; i < this->count; i++)
            this->profiles[i]->clean();
        this->resolve_dependencies();
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
