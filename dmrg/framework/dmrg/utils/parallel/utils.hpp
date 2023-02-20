/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_PARALLEL_UTILS_HPP
#define UTILS_PARALLEL_UTILS_HPP

#ifdef USE_AMBIENT
#include "utils/mkl_parallel.hpp"
#include "ambient/container/future.hpp"
#include "utils/meminfo.hpp"
#else
#include "dmrg/utils/proc_status.h"
#endif

namespace parallel {

    #ifdef USE_AMBIENT
    inline void meminfo(){
        ambient::sync();
        for(int i = 0; i < ambient::scope::size()/groups_granularity; i++){
            parallel::guard group(parallel::traits::to_iterator(i*groups_granularity), groups_granularity);
            ambient::meminfo();
        }
    }
    inline bool local(){
        return ambient::scope::local();
    }
    inline bool uniq(){
        return (ambient::master() || ambient::scope::nested());
    }
    inline void sync(){
        ambient::sync();
    }
    inline void sync_mkl_parallel(){
        ambient::sync(ambient::mkl_parallel());
    }
    inline ambient::rank_t rank(){
        return ambient::rank();
    }
    inline std::string rank_str(){
        return std::to_string(rank());
    }
    inline bool master(){
        return ambient::master();
    }
    #else
    inline void meminfo(){
        parallel::cout << "Memory usage : " << proc_status_mem() << std::endl;
    }
    inline void sync(){
    }
    inline void sync_mkl_parallel(){
    }
    inline int rank(){
        return 0;
    }
    inline std::string rank_str(){
        return "0";
    }
    inline bool master(){
        return true;
    }
    inline bool uniq(){
        return true;
    }
    inline bool local(){
        return true;
    }
    #endif

}

#endif

