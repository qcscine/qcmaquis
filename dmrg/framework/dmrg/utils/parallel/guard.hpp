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

#ifndef PARALLEL_GUARD_HPP
#define PARALLEL_GUARD_HPP

namespace parallel {

    #ifdef USE_AMBIENT
    struct guard {
        typedef ambient::actor proc_guard;
        typedef ambient::scope group_guard;
        typedef typename traits::resource_iterator iter;
        typedef typename traits::resource_iterator_nop iter_nop;

        struct serial {
            serial() : proc(ambient::actor_t::common) { }
           ~serial(){ ambient::sync(); }
        private:
            proc_guard proc;
        };

        explicit guard(const iter& r) : group(NULL) {
            proc = new proc_guard(r);
        }
        explicit guard(const iter& begin, size_t size) : proc(NULL) {
            group = new group_guard(begin, size);
        }
        explicit guard(const iter& begin, const iter& end) : proc(NULL) {
            group = new group_guard(begin, end);
        }
       ~guard(){
           if(group) delete group;
           else if(proc) delete proc;
        }
        explicit guard(iter_nop)              : proc(NULL), group(NULL) {}
        explicit guard(iter_nop, iter_nop)    : proc(NULL), group(NULL) {}
        explicit guard(iter_nop, size_t size) : proc(NULL), group(NULL) {}
    private:
        proc_guard* proc;
        group_guard* group;
    };
    #else
    struct guard {
        struct serial {};
        template<typename T> guard(const T& r){}
        template<typename T1, typename T2> guard(const T1& begin, const T2& end){}
    };
    #endif

}

#endif
