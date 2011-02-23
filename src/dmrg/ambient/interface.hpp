#ifndef AMBIENT_INTERFACE_H
#define AMBIENT_INTERFACE_H
#include <assert.h>
#include "ambient/core/operation/operation.h"
#include "ambient/core/p_profile.h"
#include "ambient/core/select.h"
#include "ambient/interface/profiles.hpp"

namespace blas{ using namespace ambient;

    template<typename T>
    void assign(T& ref, int i, int j = 0, int k = 0)
    {
        void_pt* profile = get_profile(ref);
        if(asmp.rank == UNDEFINED_RANK) return;
        workgroup* group = profile->group(i,j,k);
//        printf("%s: p%d: I've accepted group %d %d of id%d\n", asmp.get_scope()->name, asmp.rank, group->i, group->j, (*(group->profile))->id );
        group->owner = ambient::rank();
        profile->layout->update_map_entry(ambient::rank(), i, j, k); // or add_segment_entry
    }
    template<typename T>
    dim3 get_grid_dim(T& ref)
    {
        return get_profile(ref)->get_grid_dim();
    }
    template<typename T>
    dim3 get_group_dim(T& ref)
    {
        return get_profile(ref)->get_group_dim();
    }
    template<typename T>
    dim3 get_item_dim(T& ref)
    {
        return get_profile(ref)->get_item_dim();
    }
    template<typename T>
    dim3 get_group_id(T& ref)
    {
        return get_profile(ref)->get_group_id();
    }
    #include "ambient/kernels.cpp"
} 
#include "ambient/core/operation/operation.pp.hpp"
#include "ambient/interface/core.hpp"
#endif
