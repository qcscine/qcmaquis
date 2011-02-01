#include "ambient/core/operation.h"

namespace ambient{ namespace core{

    operation::operation(void(*fp)(int, int), int arg1, int arg2){
        this->operation_ptr = (void(*)())fp;
        this->arguments = (int*)malloc(sizeof(int)*2);
        this->arguments[0] = arg1;
        this->arguments[1] = arg2;
        this->prototype = &operation::prototype_duplet;
    }
    operation::operation(void(*fp)(int, int, int), int arg1, int arg2, int arg3){
        this->operation_ptr = (void(*)())fp;
        this->arguments = (int*)malloc(sizeof(int)*3);
        this->arguments[0] = arg1;
        this->arguments[1] = arg2;
        this->arguments[2] = arg3;
        this->prototype = &operation::prototype_triplet;
    }
    void operation::prototype_duplet(){ ((void(*)(int,int))this->operation_ptr)(this->arguments[0], this->arguments[1]); }
    void operation::prototype_triplet(){ ((void(*)(int,int,int))this->operation_ptr)(this->arguments[0], this->arguments[1], this->arguments[2]); }
    void operation::perform(){ (this->*prototype)(); }

} }
