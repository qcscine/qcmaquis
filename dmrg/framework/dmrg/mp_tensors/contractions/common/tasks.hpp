/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2017 Institute for Theoretical Physics, ETH Zurich
 *                    Laboratory for Physical Chemistry, ETH Zurich
 *                    Department of Chemistry and the PULSE Institute, Stanford University
 *               2017-2017 by Sebastian Keller <sebkelle@phys.ethz.ch>
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

#ifndef ENGINE_COMMON_TASKS_HPP
#define ENGINE_COMMON_TASKS_HPP

#include <vector>
#include <map>
#include <utility>

namespace contraction {
namespace common {

    namespace detail { 

        template <typename T>
        struct micro_task
        {
            typedef unsigned short IS;

            T scale;
            unsigned in_offset;
            IS b2, k, l_size, r_size, stripe, out_offset;
        };

    } // namespace detail

    template <typename T>
    struct task_compare
    {
        bool operator ()(detail::micro_task<T> const & t1, detail::micro_task<T> const & t2)
        {
            return t1.out_offset < t2.out_offset;
        }
    };

    template <class Matrix, class SymmGroup>
    struct task_capsule
    {
        typedef typename SymmGroup::charge charge;
        typedef typename Matrix::value_type value_type;
        typedef detail::micro_task<value_type> micro_task;
        typedef std::map<std::pair<charge, charge>, std::vector<micro_task>, compare_pair<std::pair<charge, charge> > > map_t;

        map_t tasks;
    };

    template <class Matrix, class SymmGroup>
    struct Schedule
    {
        typedef std::vector<contraction::common::task_capsule<Matrix, SymmGroup> > schedule_t;
    }; 
    
    template<class Matrix, class SymmGroup, class TaskCalc>
    typename Schedule<Matrix, SymmGroup>::schedule_t
    create_contraction_schedule(MPSTensor<Matrix, SymmGroup> const & initial,
                                Boundary<typename storage::constrained<Matrix>::type, SymmGroup> const & right,
                                MPOTensor<Matrix, SymmGroup> const & mpo,
                                TaskCalc task_calc)
    {
        typedef typename SymmGroup::charge charge;
        typedef typename MPOTensor<Matrix, SymmGroup>::index_type index_type;
        typedef typename Matrix::value_type value_type;
        typedef typename task_capsule<Matrix, SymmGroup>::map_t map_t;

        typename Schedule<Matrix, SymmGroup>::schedule_t contraction_schedule;

        initial.make_left_paired();

        contraction_schedule.resize(mpo.row_dim());
        contraction::common::MPSBoundaryProductIndices<Matrix, typename storage::constrained<Matrix>::type, SymmGroup>
            indices(initial.data().basis(), right, mpo);

        Index<SymmGroup> const & physical_i = initial.site_dim(),
                                 right_i = initial.col_dim();
        Index<SymmGroup> left_i = initial.row_dim(),
                         out_right_i = adjoin(physical_i) * right_i;

        common_subset(out_right_i, left_i);
        ProductBasis<SymmGroup> in_left_pb(physical_i, left_i);
        ProductBasis<SymmGroup> out_right_pb(physical_i, right_i,
                                             boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                                                 -boost::lambda::_1, boost::lambda::_2));
        index_type loop_max = mpo.row_dim();

        omp_for(index_type b1, parallel::range<index_type>(0,loop_max), {
            task_capsule<Matrix, SymmGroup> tasks_cap;
            task_calc(b1, indices, mpo, initial.data().basis(), left_i, out_right_i, in_left_pb, out_right_pb, tasks_cap);

            map_t & tasks = tasks_cap.tasks;
            for (typename map_t::iterator it = tasks.begin(); it != tasks.end(); ++it)
                std::sort((it->second).begin(), (it->second).end(), contraction::common::task_compare<value_type>());

            contraction_schedule[b1] = tasks_cap;
        });

        size_t sz = 0;
        for (int b1 = 0; b1 < loop_max; ++b1)
        {
            map_t & tasks = contraction_schedule[b1].tasks;
            for (typename map_t::iterator it = tasks.begin(); it != tasks.end(); ++it)
                sz += (it->second).size() * sizeof(contraction::common::detail::micro_task<value_type>);
        }
        maquis::cout << "Schedule size: " << sz / 1024 << std::endl;

        return contraction_schedule;
    }


} // namespace common
} // namespace contraction

#endif
