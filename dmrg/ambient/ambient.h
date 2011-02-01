#ifndef AMBIENT_H
#define AMBIENT_H

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <queue>

#include "ambient/groups/group.h"
#include "ambient/groups/multirank.h"
#include "ambient/auxiliary.h"
#include "ambient/core/operation.h"

namespace ambient
{
    class scheduler
    {
    private: 
        scheduler();                             // constructor is private
        scheduler(scheduler const&);             // copy constructor is private
        scheduler& operator=(scheduler const&);  // assignment operator is private
    public:
        static scheduler& instance();

    public:
        scheduler & operator>>(dim3 dim_distr);
        scheduler & operator,(dim3 dim);
        bool is_ambient_master(); 
        void init(MPI_Comm comm = NULL);
        void regression_test();
        void finalize();
        dim3 group_dim();
        dim3 item_dim();

        void push(ambient::core::operation* logistics, ambient::core::operation* computing);
        void evaluate_op_stack(); 
        void perform_op_stack(); 
    private:
        MPI_Comm comm;
        int size;
        groups::group* ambient;
        groups::multirank& rank;

        enum { AMBIENT_MASTER,
               GROUP_MASTER,
               GROUP_SLAVE } mode;

        dim3 dim_distr;   // work-item size of distribution blocks
        dim3 dim_group;   // work-item size of cpu streaming multiprocessor workload fractions
        dim3 dim_item;    // size of work-item (i.e. 128) 
        dim3 dim_gpu;     // work-item size of gpgpu smp workload fractions

        std::queue<ambient::core::operation*> logistics_stack;
        std::queue<ambient::core::operation*> computing_stack;
    };

    scheduler& operator>>(scheduler* instance, dim3 dim_distr);
    size_t get_bound();
    scheduler& instance();
}

#endif
