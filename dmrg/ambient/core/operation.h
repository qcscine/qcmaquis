#ifndef AMBIENT_CORE_OPERATION_H
#define AMBIENT_CORE_OPERATION_H
#include <stdlib.h>
#include "ambient/core/p_profile.h"

namespace ambient{ namespace core{

    typedef p_profile   void_pt;
    typedef p_profile_s void_spt;

    class operation{
    public:
        operation(void(*fp)(void_pt&, void_pt&, void_spt&), void_pt* arg1, void_pt* arg2, void_spt* arg3);
        operation(void(*fp)(const void_pt&, void_pt&, void_spt&), void_pt* arg1, void_pt* arg2, void_spt* arg3);
        void prototype_triplet();
        void prototype_triplet_const();
        void perform(); // executes operation
        void performx(); // executes operation (clean way)
        void set_ids();
        void set_scope(groups::group* scope);
        groups::group* get_scope();
        void(operation::*prototype)();
        void(*operation_ptr)();
        groups::group* scope;
        void_pt** arguments;
        void_spt* pin;
        size_t arg_count;
    };

    class out_of_scope_e{};
    class interrupt_e{};

} }

#endif
