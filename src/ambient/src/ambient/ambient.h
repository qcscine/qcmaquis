#ifndef AMBIENT_H
#define AMBIENT_H

#define _AMBIENT

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <queue>
#include <string>
#include <assert.h>
#include <limits> // C - for the validation test

#include "ambient/traits.h"
#include "ambient/groups/group.h"
#include "ambient/groups/multirank.h"
#include "ambient/auxiliary.h"
#include "ambient/core/scope_context.h"
#include "ambient/core/operation/operation.h"
#include "ambient/core/workgroup_context.h"
#include "ambient/model.h"
#include "ambient/core/layout.h"
#include "ambient/core/select.h"
#include "ambient/core/auxiliary.h"
#include <boost/shared_ptr.hpp>

class Timer;

#define ALL -1
#define UNDEFINED_RANK MPI_UNDEFINED
#define ID_TYPE unsigned long long int

#define __ambient_wo_begin__ ambient::access.write_only_start_mark();
#define __ambient_wo_end__   ambient::access.write_only_stop_mark();


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
        scheduler & operator>>(dim2 mem_dim);
        scheduler & operator,(dim2 dim);
        void init();
        void finalize();
        dim2 get_mem_dim();
        dim2 get_item_dim();
        dim2 get_work_dim();
        dim2 get_gpu_dim();

        void push(core::operation* logistics, core::operation* computing);

        inline void playout();  // perform actual operations
        inline void bailout();  // perform actual operations
        inline void bailin();  // perform actual operations
        inline void spin();
        inline void spin_loop();
        inline void world_loop();
        inline void world_spin();
        inline bool occupied();
        inline bool blank();

        int size;
        block_packet_t<typename traits::type>* default_data_packet_t;
        groups::group* ambient;
    private:
        dim2 work_dim;   // work-item size of distribution blocks
        dim2 mem_dim;     // work-item size of cpu streaming multiprocessor workload fractions
        dim2 item_dim;    // size of work-item (i.e. 128) 
        dim2 gpu_dim;     // work-item size of gpgpu smp workload fractions
        bool stirring;    // playout is active

        one_touch_stack< std::pair<core::operation*,core::operation*> > stack;
        one_touch_stack< groups::packet_manager* > router; // packet_manager router
        workgroup_context context;
    };

    scheduler& operator>>(scheduler* instance, dim2 mem_dim);
    void init();
    void finalize();
    void playout();
    void bailout();
    void bailin();
    void spin();
    void spin_loop();
    void world_loop();
    void world_spin();
    int size();
    bool is_master();
    bool is_group_master();
    bool occupied();
    bool blank();
    groups::group* world();

    extern scope_context& scope;
    extern scheduler& layout;
    extern scheduler& engine;
    extern groups::multirank& rank;
    extern i_model&  model;
    extern i_controller& controller;
    extern i_channel& channel;
    extern groups::comm_map& scope_map;
    extern access_marker& access;
    extern Timer timer;
}

#endif
