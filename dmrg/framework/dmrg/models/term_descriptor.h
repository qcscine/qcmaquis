/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 *               2021 by Alberto Baiardi <abaiardi@ethz.ch>
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

#ifndef MODELS_TERM_DESCRIPTOR_H
#define MODELS_TERM_DESCRIPTOR_H

#include <utility>
#include <vector>

namespace detail {
    
struct pos_tag_lt {
    using value_type = std::pair<int, unsigned int>;
    inline bool operator() (value_type const& lhs, value_type const& rhs)
    {
        return (lhs.first < rhs.first);
    }
};
    
}

/**
 * @brief Term descriptor class
 * 
 * The term descriptor class represents an entry of an operators expressed in second-quantization.
 * The underlying data type is a vector of pairs (int, int). The first element of the pair is
 * the site on which the operator acts, and the second element is the tag of the corresponding elementary
 * operator.
 * 
 * @tparam T Scalar type underlying the operator.
 */
template <typename T>
class term_descriptor : public std::vector<std::pair<int, unsigned int> > {
public:
    // Types definition
    using pos_type = int;
    using tag_type = unsigned int;
    using value_type = std::pair<pos_type, tag_type>;
    using base = std::vector<value_type>;
    using size_type = base::size_type;
    using iterator = typename base::iterator;
    using const_iterator = typename base::const_iterator;
    
    /** @brief Class constructor */
    term_descriptor() : base(), coeff(1.), is_fermionic(false) { }
    
    /** @brief Getter for the position */
    pos_type position(size_type i) const { return this->operator[](i).first; }

    /** @brief Getter for the operator */
    tag_type operator_tag(size_type i) const { return this->operator[](i).second; }
    
    /** @brief Sorting of the operators */
    // TODO: check and fix for fermions
    void canonical_order() { std::sort(begin(), end(), detail::pos_tag_lt()); }
    
    /** @brief Comparison operator */
    bool operator< (term_descriptor const & rhs) const
    {
        if (this->size() == 0 && rhs.size() == 0)
            return false;
        if (this->size() == 0)
            return true;
        if (rhs.size()   == 0)
            return false;
        if (this->position(0) == rhs.position(0))
            return this->size() > rhs.size();
        return this->position(0) < rhs.position(0);
    }

    // Class members
    /** Scaling factor for the operator term */
    T coeff;
    /** True for fermionic operators */
    bool is_fermionic;
};

// ostream
template<typename T>
std::ostream & operator<< (std::ostream & os, term_descriptor<T> const& term)
{
    os << "coeff: " << term.coeff << std::endl;
    os << "operators:";
    for (int i=0; i<term.size(); ++i)
        os << " {"  << term.position(i) << "," << term.operator_tag(i) << "}";
    os << std::endl;
    return os;
}


#endif
