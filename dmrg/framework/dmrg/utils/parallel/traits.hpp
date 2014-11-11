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

#ifndef PARALLEL_TRAITS_HPP
#define PARALLEL_TRAITS_HPP

#include <vector>
#include <cassert>

namespace parallel {

    #ifdef USE_AMBIENT
    struct traits {
        typedef int resource_iterator_nop;
        typedef typename ambient::scope::const_iterator resource_iterator;
        static resource_iterator balance(int k, int max){ return ambient::scope::balance(k, max); }
        static resource_iterator permute(int k, const std::vector<int>& s, size_t gran){ if(k >= s.size()) return to_iterator(k % size()); return ambient::scope::permute(k, s, gran); }
        static size_t size(){ return ambient::scope::size(); }
        static size_t distance(const resource_iterator& it){ return (it - ambient::scope::begin()); }
        template<typename T>
        static resource_iterator to_iterator(T offset){ assert(offset < size()); return (ambient::scope::begin() + offset); }
        template<typename T>
        static int placement(const T& obj){ return ambient::get_owner(obj); }
    };
    #else
    struct traits {
        typedef int resource_iterator_nop;
        typedef resource_iterator_nop resource_iterator;
        static resource_iterator balance(int k, int max){ return resource_iterator(); }
        static resource_iterator permute(int k, const std::vector<int>& s, size_t){ return resource_iterator(); }
        static size_t size(){ return 1; }
        template<typename T>
        static resource_iterator to_iterator(T offset){ return resource_iterator(); }
        template<typename T>
        static int placement(const T& obj){ return 0; }
    };
    #endif

}

#endif
