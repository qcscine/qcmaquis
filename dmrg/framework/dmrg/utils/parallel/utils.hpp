/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2013 by Alexandr Kosenkov <alex.kosenkov@gmail.com>
 *                         by Michele Dolfi <dolfim@phys.ethz.ch>
 * 
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 * 
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not, the license is also
 * available from http://alps.comp-phys.org/.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
 * DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

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

