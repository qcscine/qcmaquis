/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
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

#ifndef PARALLEL_FOR_HPP
#define PARALLEL_FOR_HPP

#ifdef AMBIENT
    typedef ambient::scope<ambient::scope_t::single> locale;
    typedef ambient::scope<ambient::scope_t::shared> locale_shared;
    #define parallel_for(constraint, ...) constraint; for(__VA_ARGS__)
    #define semi_parallel_for(constraint, ...) constraint; for(__VA_ARGS__)
#elif defined(MAQUIS_OPENMP)
    typedef std::size_t locale;
    typedef std::size_t locale_shared;
    #define parallel_pragma(a) _Pragma( #a )
    #define parallel_for(constraint, ...) parallel_pragma(omp parallel for schedule(dynamic, 1)) for(__VA_ARGS__)
    #define semi_parallel_for(constraint, ...) for(__VA_ARGS__)
#else
    typedef std::size_t locale;
    typedef std::size_t locale_shared;
    #define parallel_for(constraint, ...) for(__VA_ARGS__)
    #define semi_parallel_for(constraint, ...) for(__VA_ARGS__)
#endif

#endif
