/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2014-2014 by Sebastian Keller <sebkelle@phys.ethz.ch>
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

#ifndef SPIN_DESCRIPTOR_H
#define SPIN_DESCRIPTOR_H

#include "dmrg/block_matrix/symmetry/symmetry_traits.h"

template <class SymmType>
class SpinDescriptor
{
public:
    typedef int spin_t;

    void clear() { }
    int get() const { return 0; }

    SpinDescriptor() {}
    SpinDescriptor(spin_t twoS_, spin_t in, spin_t out) {}

    bool operator==(SpinDescriptor const & rhs) const { return true; }
    bool operator!=(SpinDescriptor const & rhs) const { return false; }

    template <class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {}
};

template <class SymmType> inline
SpinDescriptor<SymmType> couple(SpinDescriptor<SymmType> a, SpinDescriptor<SymmType> b)
{
    return SpinDescriptor<SymmType>();
}

template <class SymmType> inline
std::ostream & operator << (std::ostream & os, SpinDescriptor<SymmType> rhs) { return os; }

template <>
class SpinDescriptor<symm_traits::SU2Tag>
{
public:
    typedef int spin_t;

    SpinDescriptor() : twoS(0), diff_(0) {}
    SpinDescriptor(spin_t twoS_, spin_t in, spin_t out) : twoS(twoS_), diff_(out-in) {}

    SpinDescriptor & operator += (SpinDescriptor rhs)
    {
        // apply action of operator rhs 
        twoS += rhs.action();
        return *this;
    }

    spin_t get() const { return twoS; }
    spin_t action() const { return diff_; }

    void clear()
    {
        twoS = 0; diff_ = 0;
    }

    bool operator==(SpinDescriptor const & rhs) const { return (twoS == rhs.twoS && diff_ == rhs.diff_); }
    bool operator!=(SpinDescriptor const & rhs) const { return !(*this==rhs); }
    friend SpinDescriptor operator-(SpinDescriptor);
    friend std::ostream & operator<<(std::ostream & os, SpinDescriptor);

    template <class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & twoS & diff_;
    }

private:
    spin_t twoS;
    spin_t diff_;
};

// Attention: not symmetric
inline SpinDescriptor<symm_traits::SU2Tag> couple(SpinDescriptor<symm_traits::SU2Tag> a, SpinDescriptor<symm_traits::SU2Tag> b)
{
    return a += b;
}

inline SpinDescriptor<symm_traits::SU2Tag> operator-(SpinDescriptor<symm_traits::SU2Tag> rhs)
{
    rhs.diff_ = -rhs.diff_;
    return rhs;
}

inline std::ostream & operator<<(std::ostream & os, SpinDescriptor<symm_traits::SU2Tag> rhs)
{
    os << "Spin: " << rhs.get() << ", delta: " << rhs.action() << std::endl;
    return os;
}


// This function only works if we are dealing with spin 1/2 particles
// In this case we can map the product spin 1 to 1-1×11 and spin 0 to 11×1-1

template <class SymmGroup>
typename SymmGroup::subcharge productSpin(typename SymmGroup::charge a, typename SymmGroup::charge b)
{
    typename SymmGroup::subcharge spin_a = SymmGroup::spin(a);
    typename SymmGroup::subcharge spin_b = SymmGroup::spin(b);

    if (spin_a == -1 && spin_b == 1)
        return 2;
    else
        return std::abs(spin_a + spin_b);
}

#endif
