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

#ifndef PARALLEL_SCHEDULER_BALANCED_HPP
#define PARALLEL_SCHEDULER_BALANCED_HPP

namespace parallel {

    class scheduler_balanced {
    public:
        typedef traits::resource_iterator resource_iterator;

        template<class Matrix>
        scheduler_balanced(const Matrix& m) : max_k(m.n_blocks()) {}
        scheduler_balanced(size_t max) : max_k(max) {}
        resource_iterator operator()(int k) const {
            return traits::balance(k,max_k);
        }
        bool propagate() const {
            return false;
        }
    protected:
        int max_k;
    };

}

#endif
