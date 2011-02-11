#ifndef AMBIENT_H
#define AMBIENT_H

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <queue>
#include <string>

#include "ambient/groups/group.h"
#include "ambient/groups/multirank.h"
#include "ambient/auxiliary.h"
#include "ambient/core/smp.h"
#include "ambient/core/operation.h"
#include "ambient/core/coherency.h"

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
        bool is_master() const; 
        bool is_ambient_master() const; 
        void init(MPI_Comm comm = NULL);
        void regression_test();
        void finalize();
        dim3 group_dim();
        dim3 item_dim();

        void push(ambient::core::operation* logistics, ambient::core::operation* computing);
        void playout();  // perform actual operations
        int size;
    private:
        MPI_Comm comm;
        groups::group* ambient;

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
    void playout();
    int size();

    extern smp& asmp;
    extern scheduler& layout;
    extern scheduler& engine;
    extern groups::multirank& rank;
    extern core::coherency_table& coherency;
    extern hash_map<core::coherency_table> void_pt_map;
}

#endif
