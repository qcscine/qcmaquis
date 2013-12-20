/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef MEASUREMENT_H
#define MEASUREMENT_H

#include "dmrg/mp_tensors/reduced_mps.h"
#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"
#include "dmrg/block_matrix/symmetry.h"

#include "dmrg/models/lattice.h"

#include <vector>
#include <string>
#include <sstream>
#include <iostream>

#include <stdexcept>

#include <boost/regex.hpp>
#include <boost/static_assert.hpp>

/// TODO: 1) move data to new object measurement_result.
///       2) store only measurement description, i.e. no matrix, and pass model+lattice in evaluate method.

template<class Matrix, class SymmGroup>
class measurement {
public:
    typedef typename Matrix::value_type value_type;
    
    measurement(std::string const& n="") : cast_to_real(true), is_super_meas(false), name_(n) { }
    virtual ~measurement() { }
    
    virtual void evaluate(MPS<Matrix, SymmGroup> const&, boost::optional<reduced_mps<Matrix, SymmGroup> const&> = boost::none) =0;
    template <class Archive>
    void save(Archive &) const;
    virtual void print(std::ostream& os) const;
    
    std::string const& name() const { return name_; }
    
    void set_super_meas(Index<SymmGroup> const& phys_psi_);
    
    measurement* clone() const { return do_clone(); }
    
protected:
    virtual measurement* do_clone() const =0;
    
    bool cast_to_real, is_super_meas;
    std::vector<std::string> labels;
    value_type result;
    std::vector<value_type> vector_results;
    
    Index<SymmGroup> phys_psi;
    
private:
    std::string name_;
};

template<class Matrix, class SymmGroup>
inline measurement<Matrix, SymmGroup>* new_clone( const measurement<Matrix, SymmGroup>& m )
{
    return m.clone();
}

template<class Matrix, class SymmGroup>
template <class Archive>
void measurement<Matrix, SymmGroup>::save(Archive & ar) const
{
    if (labels.size() > 0) {
        if (cast_to_real) {
            std::vector<std::vector<double> > tmp(1, maquis::real(vector_results));
            ar[storage::encode(name()) + std::string("/mean/value")] << tmp;
        } else {
            ar[storage::encode(name()) + std::string("/mean/value")] << std::vector<std::vector<value_type> >(1, vector_results);
        }
        ar[storage::encode(name()) + std::string("/labels")] << labels;
    } else {
        if (cast_to_real) {
            std::vector<double> tmp(1, maquis::real(result));
            ar[storage::encode(name()) + std::string("/mean/value")] << tmp;
        } else {
            ar[storage::encode(name()) + std::string("/mean/value")] << std::vector<value_type>(1, result);
        }
    }
}

template<class Matrix, class SymmGroup>
void measurement<Matrix, SymmGroup>::set_super_meas(Index<SymmGroup> const& phys_psi_)
{
    phys_psi = phys_psi_;
    is_super_meas = true;
}

template<class Matrix, class SymmGroup>
void measurement<Matrix, SymmGroup>::print(std::ostream& os) const
{
    os << "MEASURE[" << name_ << "]";
}

template<class Matrix, class SymmGroup>
std::ostream& operator<<(std::ostream& os, measurement<Matrix, SymmGroup> const& m)
{
    m.print(os);
}


/// UTILITIES

template<class Matrix, class SymmGroup>
bool is_hermitian_meas(std::vector<block_matrix<Matrix, SymmGroup> > const & ops)
{
    return all_true(ops.begin(), ops.end(),
                    boost::bind(static_cast<bool (*)(block_matrix<Matrix, SymmGroup> const&)>(&is_hermitian), _1));
}

template<class Matrix, class SymmGroup>
bool is_hermitian_meas(std::vector<std::pair<std::vector<block_matrix<Matrix, SymmGroup> >, bool> > const & ops)
{
    bool is_herm = true;
    for (int i=0; i<ops.size() && is_herm; ++i)
        is_herm = is_hermitian_meas(ops[i].first);
    return is_herm;
}

inline std::vector<std::string> label_strings (const Lattice& lat, const std::vector<std::vector<std::size_t> >& labels)
{
    std::vector<std::string> ret;
    ret.reserve(labels.size());
    for (std::vector<std::vector<std::size_t> >::const_iterator it = labels.begin();
         it != labels.end(); ++it)
    {
        std::ostringstream oss;
        for (std::vector<std::size_t>::const_iterator it2 = it->begin(); it2 != it->end(); ++it2) {
            oss << lat.get_prop<std::string>("label", *it2);
            if (it2 + 1 != it->end())
                oss << " -- ";
        }
        ret.push_back(oss.str());
    }
    return ret;
}



#endif
