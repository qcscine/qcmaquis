#ifndef AMBIENT_H
#define AMBIENT_H

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "dim3.h"
#include "ambient/groups/group.h"
#include "ambient/groups/multirank.h"
#include "ambient/interface/p_profile.h"

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

    private:
        MPI_Comm comm;
        int size;
        groups::group* ambient;
        groups::multirank& rank;

        enum { AMBIENT_MASTER,
               GROUP_MASTER,
               GROUP_SLAVE } mode;

        dim3 dim_distr; // work-item size of distribution blocks
        dim3 dim_cpu;   // work-item size of cpu streaming multiprocessor workload fractions
        dim3 dim_gpu;   // work-item size of gpgpu smp workload fractions
    };

    scheduler& operator>>(scheduler* instance, dim3 dim_distr);
    size_t get_bound();
}

#endif
