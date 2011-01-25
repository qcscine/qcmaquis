#ifndef AMBIENT_P_PROFILE_H
#define AMBIENT_P_PROFILE_H
#include "ambient/dim3.h"
#include "ambient/interface/p_action.h"
#include <utility>

namespace ambient {

    class p_profile {
    public:
        template <typename T>
        p_profile(const T* ptr);
        const void* ptr;
        const char* type;
        p_action* action;
        bool proxy;
        dim3 size;
        dim3 block_size;
        int** owners;
    };

}
#endif
