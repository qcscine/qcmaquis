#ifndef AMBIENT_P_PROFILE_H
#define AMBIENT_P_PROFILE_H
#include "ambient/dim3.h"
#include "ambient/interface/p_action.h"
#include <utility>
#include <list>
#include <vector>

namespace ambient {

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
        
        dim3 dim_block;
        dim3 dim_group;
        int** owners;

        p_action* action; // do I need this -- apparently not

        std::vector<workgroup*> skeleton;
        workgroup* group(int i, int j = 0, int k = 0);
    };

}
#endif
