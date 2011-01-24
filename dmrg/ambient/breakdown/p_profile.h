#ifndef AMBIENT_P_PROFILE_H
#define AMBIENT_P_PROFILE_H
#include "ambient/dim3.h"
#include <utility>

namespace ambient { namespace breakdown {

    class p_profile {
    public:
        p_profile(const void* ptr, const char* type, dim3 size, dim3 block_size, int** owners)
        : ptr(ptr), type(type), size(size), block_size(block_size), owners(owners) { };
        p_profile(const void* ptr, const char* type) : ptr(ptr), type(type) { };
        const void* ptr;
        const char* type;
        dim3 size;
        dim3 block_size;
        int** owners;
    };

} }
#endif
