#ifndef AMBIENT_P_PROFILE_H
#define AMBIENT_P_PROFILE_H
#include "ambient/dim3.h"
#include <utility>
#include <list>
#include <vector>

namespace ambient {

    class p_profile;

    class workgroup {
    public:
        workgroup(p_profile* p, int i, int j);
        void* item(int i, int j = 0, int k = 0);
        p_profile* profile;
        int owner;
        int i, j;
    };

    class p_profile {
    public:
        template <typename T>
        p_profile(const T* ptr);

        void* scope; // includes ambient boundle
        void* data;  // pointer to the actual data
        size_t lda;  // process individual lda

        const char* type;
        bool proxy;
        int** owners;

        std::vector<workgroup*> skeleton;
        workgroup* group(int i, int j = 0, int k = 0);
    };

}
#endif
