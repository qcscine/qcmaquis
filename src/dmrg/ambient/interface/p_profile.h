#ifndef AMBIENT_P_PROFILE_H
#define AMBIENT_P_PROFILE_H
#include "ambient/interface/workgroup.h"
#include "ambient/auxiliary.h"
#include <utility>
#include <list>
#include <vector>

namespace ambient {

    class workgroup;

    class p_profile {
    public:
        int id;
        template <typename T>
        p_profile(const T* ptr);
        p_profile & operator>>(dim3 dim_distr);
        p_profile & operator,(dim3 dim);
        void regroup();

        p_profile* profile; // pointer to this profile (this on init - can be changed in proxy objects)
        p_profile* dereference(); // finds out if the profile pointer is up to date
        void* scope; // includes ambient boundle
        void* data;  // pointer to the actual data
        size_t lda;  // process individual lda

        const char* type;
        bool proxy;
        
        dim3 dim;
        int** owners;

        std::vector<workgroup*> skeleton;
        workgroup* group(int i, int j = 0, int k = 0);

        dim3 grid_dim();
        dim3 group_dim();
        dim3 item_dim();
// parameters can be set specifically for the profile
        bool specific; 
        dim3 dim_distr;   // work-item size of distribution blocks
        dim3 dim_group;   // work-item size of cpu streaming multiprocessor workload fractions
        dim3 dim_item;    // size of work-item (i.e. 128) 
        dim3 dim_gpu;     // work-item size of gpgpu smp workload fractions
    };

    p_profile& operator>>(p_profile* instance, dim3 dim_distr);

}
#endif
