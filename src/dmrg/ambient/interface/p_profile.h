#ifndef AMBIENT_P_PROFILE_H
#define AMBIENT_P_PROFILE_H
#include "ambient/ambient.h"
#include "ambient/interface/workgroup.h"
#include <utility>
#include <list>
#include <vector>

namespace ambient {

    class workgroup;

    class p_profile {
    public:
        template <typename T>
        p_profile(const T* ptr);

        void* scope; // includes ambient boundle
        void* data;  // pointer to the actual data
        size_t lda;  // process individual lda

        const char* type;
        bool proxy;
        
        dim3 dim;
        int** owners;

        std::vector<workgroup*> skeleton;
        workgroup* group(int i, int j = 0, int k = 0);
    };

}
#endif
