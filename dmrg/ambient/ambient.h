#ifndef AMBIENT_H
#define AMBIENT_H

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <queue>
#include <string>
#include <assert.h>

#include "ambient/groups/group.h"
#include "ambient/groups/multirank.h"
#include "ambient/auxiliary.h"
#include "ambient/core/smp.h"
#include "ambient/core/operation/operation.h"
#include "ambient/core/layout.h"
#include "ambient/core/select.h"

#define ALL -1
#define UNDEFINED_RANK MPI_UNDEFINED
#define ID_TYPE unsigned long long int

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
        void init(MPI_Comm comm = NULL);
        void regression_test();
        void finalize();
        dim3 group_dim();
        dim3 item_dim();

        void push(core::operation* logistics, core::operation* computing);
        void playout();  // perform actual operations
        int size;
    private:
        MPI_Comm comm;
        groups::group* ambient;

        dim3 dim_distr;   // work-item size of distribution blocks
        dim3 dim_group;   // work-item size of cpu streaming multiprocessor workload fractions
        dim3 dim_item;    // size of work-item (i.e. 128) 
        dim3 dim_gpu;     // work-item size of gpgpu smp workload fractions

        std::list<core::operation*> logistics_stack;
        std::list< std::pair<core::operation*,core::operation*> > computing_stack;
    };

    scheduler& operator>>(scheduler* instance, dim3 dim_distr);
    size_t get_bound();
    size_t get_block_bound();
    void playout();
    int size();
    bool is_master();

    extern smp& asmp;
    extern scheduler& layout;
    extern scheduler& engine;
    extern groups::multirank& rank;
    extern hash_map& p_profile_map;
}

#endif
