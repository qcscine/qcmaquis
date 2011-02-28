#ifndef AMBIENT_CORE_SMP_H
#define AMBIENT_CORE_SMP_H

#include "ambient/core/operation/operation.h"
#include "ambient/core/p_profile.h"
#include "ambient/core/workgroup.h"
#include "ambient/core/layout.h"

namespace ambient {

    class smp { // scalable multiprocessor
    private:
        smp();
        smp(smp const&);
        smp& operator=(smp const&);
    public:
        static smp& instance();

    public:
        smp& operator()(const int rank);
        void set_group(groups::group* group);
        void set_group(const char* group);
        groups::group* get_group();
        int get_rank();
        int get_size();
        void set_op(core::operation* op);
        core::operation* get_op();
        bool involved();
        bool master(); // if the invoking process is master of group
    private:
        int rank;
        int size;
        groups::group* group;
        core::operation* op;
        core::layout_table* assignment;
        std::list<workgroup*>  sendlist;
        std::list<workgroup*>  recvlist;
    };

}
#endif
