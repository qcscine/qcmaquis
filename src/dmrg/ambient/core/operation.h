#ifndef AMBIENT_CORE_OPERATION_H
#define AMBIENT_CORE_OPERATION_H

#include <stdlib.h>

namespace ambient{ namespace core{

    class operation{
    public:
        operation(void(*fp)(int, int), int arg1, int arg2);
        operation(void(*fp)(int, int, int), int arg1, int arg2, int arg3);
        void prototype_duplet();
        void prototype_triplet();
        void perform(); // executes operation
        void(operation::*prototype)();
        void(*operation_ptr)();
        int* arguments;
    };

} }

#endif
