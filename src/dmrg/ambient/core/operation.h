#ifndef AMBIENT_CORE_OPERATION_H
#define AMBIENT_CORE_OPERATION_H

#include <stdlib.h>
#include "ambient/interface/p_profile.h"

namespace ambient{ namespace core{

    class operation{
    public:
        operation(void(*fp)(p_profile*, p_profile*), p_profile* arg1, p_profile* arg2);
        operation(void(*fp)(p_profile*, p_profile*, p_profile*), p_profile* arg1, p_profile* arg2, p_profile* arg3);
        void prototype_duplet();
        void prototype_triplet();
        void perform(); // executes operation
        void(operation::*prototype)();
        void(*operation_ptr)();
        p_profile** arguments;
    };

} }

#endif
