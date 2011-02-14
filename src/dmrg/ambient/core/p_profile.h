#ifndef AMBIENT_CORE_P_PROFILE_H
#define AMBIENT_CORE_P_PROFILE_H
#include "ambient/core/workgroup.h"
#include "ambient/auxiliary.h"
#include <utility>
#include <list>
#include <vector>

namespace ambient {

    class workgroup;
    class p_profile;

    class p_profile_s {
    protected:
        p_profile_s();
    public:
        unsigned int* group_id;
        unsigned int id;
        p_profile & operator>>(dim3 dim_distr);
        p_profile & operator,(dim3 dim);
        void regroup();
        void set_id(std::pair<unsigned int*,size_t> group_id);

        p_profile* profile; // pointer to this profile (this on init - can be changed in proxy objects)
        p_profile* dereference(); // finds out if the profile pointer is up to date
        void* scope; // includes ambient boundle
        void* data;  // pointer to the actual data
        size_t lda;  // process individual lda
        core::layout_table* layout;

        size_t reserved_x;
        size_t reserved_y;

        const char* type;
        bool proxy;
        
        dim3 dim;
        int** owners;

        std::vector< std::vector<workgroup*> > skeleton;

        void*(*init_fp)(workgroup* grp);
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

    class p_profile: public p_profile_s 
    { public: template <typename T> p_profile(const T* ptr, void*(*init_fp)(workgroup* grp) = NULL); };

    p_profile& operator>>(p_profile* instance, dim3 dim_distr);

}
#endif
