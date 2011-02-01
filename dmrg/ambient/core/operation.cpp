#include "ambient/core/operation.h"

namespace ambient{ namespace core{

    operation::operation(void(*fp)(p_profile*, p_profile*), p_profile* arg1, p_profile* arg2){
        this->operation_ptr = (void(*)())fp;
        this->arguments = (p_profile**)malloc(sizeof(p_profile*)*2);
        this->arguments[0] = arg1;
        this->arguments[1] = arg2;
        this->prototype = &operation::prototype_duplet;
    }
    operation::operation(void(*fp)(p_profile*, p_profile*, p_profile*), p_profile* arg1, p_profile* arg2, p_profile* arg3){
        this->operation_ptr = (void(*)())fp;
        this->arguments = (p_profile**)malloc(sizeof(p_profile*)*3);
        this->arguments[0] = arg1;
        this->arguments[1] = arg2;
        this->arguments[2] = arg3;
        this->prototype = &operation::prototype_triplet;
    }
    void operation::prototype_duplet(){ ((void(*)(p_profile*,p_profile*))this->operation_ptr)(this->arguments[0], this->arguments[1]); }
    void operation::prototype_triplet(){ ((void(*)(p_profile*,p_profile*,p_profile*))this->operation_ptr)(this->arguments[0], this->arguments[1], this->arguments[2]); }
    void operation::perform(){ (this->*prototype)(); }

} }
