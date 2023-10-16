/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef PARALLEL_GUARD_HPP
#define PARALLEL_GUARD_HPP

namespace parallel {

    #ifdef USE_AMBIENT
    struct guard {
        typedef ambient::actor proc_guard;
        typedef ambient::actor_common proc_guard_common;
        typedef ambient::scope group_guard;
        typedef typename traits::resource_iterator iter;
        typedef typename traits::resource_iterator_nop iter_nop;

        struct serial {
           ~serial(){ ambient::sync(); }
        private:
            proc_guard_common proc;
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
