/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#include <cmath>
#include <iterator>
#include <iostream>

using std::cerr;
using std::cout;
using std::endl;

#include "dmrg/block_matrix/detail/alps.hpp"

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/mp_tensors/mpstensor.h"


typedef U1 SymmGroup;
typedef alps::numeric::matrix<double> matrix;


int main() {
    
    int Nrep = 10;
    int M = 50;
    
    // Bosons with Nmax=2
    Index<SymmGroup> phys;
    phys.insert(std::make_pair(0, 1));
    phys.insert(std::make_pair(1, 1));
    phys.insert(std::make_pair(2, 1));
    
    for (int i=0; i<Nrep; ++i) {
        // create random left_i
        Index<SymmGroup> left_i;
        size_t max_left = long(dmrg_random::uniform()*10)+1;
        for (size_t k = 0; k<max_left; ++k)
            left_i.insert(std::make_pair(k, long(dmrg_random::uniform()*M)+1));
        
        // create consistent random right_i
        Index<SymmGroup> right_i = phys * left_i;
        for (long k=right_i.size()-1; k>=0; --k) {
            if (dmrg_random::uniform() < 0.2)
                right_i.erase(right_i.begin()+k);
            else
                right_i[k].second = long(dmrg_random::uniform()*M)+1);
        }
        
        MPSTensor<matrix, SymmGroup> m1(phys, left_i, right_i), m2(phys, left_i, right_i);
        
        double norm = m1.scalar_norm();
        double overlap = m1.scalar_overlap(m2);
        
    }

}
