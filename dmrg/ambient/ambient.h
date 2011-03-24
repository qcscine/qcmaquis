#ifndef AMBIENT_H
#define AMBIENT_H

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <queue>
#include <string>
#include <assert.h>

#include "ambient/groups/group.h"
#include "ambient/groups/multirank.h"
#include "ambient/auxiliary.h"
#include "ambient/core/scope_proxy.h"
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
        scheduler & operator>>(dim3 distr_dim);
        scheduler & operator,(dim3 dim);
        void init(MPI_Comm comm = NULL);
        void regression_test();
        void finalize();
        dim3 get_group_dim();
        dim3 get_item_dim();
        dim3 get_distr_dim();
        dim3 get_gpu_dim();

        void push(core::operation* logistics, core::operation* computing);
        void playout();  // perform actual operations
        void spin();
        void spin_loop();
        int size;
        block_packet_t* default_data_packet_t;
        groups::group* ambient;
    private:
        MPI_Comm comm;

        dim3 distr_dim;   // work-item size of distribution blocks
        dim3 group_dim;   // work-item size of cpu streaming multiprocessor workload fractions
        dim3 item_dim;    // size of work-item (i.e. 128) 
        dim3 gpu_dim;     // work-item size of gpgpu smp workload fractions

        one_touch_stack< std::pair<core::operation*,core::operation*> > stack;
        one_touch_stack< groups::packet_manager* > router; // packet_manager router
    };

    scheduler& operator>>(scheduler* instance, dim3 distr_dim);
    void init(MPI_Comm comm = NULL);
    void finalize();
    void playout();
    void spin();
    void spin_loop();
    int size();
    bool is_master();
    groups::group* world();

    extern scope_proxy& scope;
    extern scheduler& layout;
    extern scheduler& engine;
    extern groups::multirank& rank;
    extern hash_map& p_profile_map;
}

#endif
