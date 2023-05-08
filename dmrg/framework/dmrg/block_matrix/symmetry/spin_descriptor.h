/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

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
    int action() const { return 0; }

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

    spin_t get() const { return twoS; }
    spin_t action() const { return diff_; }

    void clear()
    {
        twoS = 0; diff_ = 0;
    }

    SpinDescriptor & couple (SpinDescriptor rhs)
    {
        // apply action of operator rhs 
        twoS += rhs.action();
        //diff_ += rhs.action();
        return *this;
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
    return a.couple(b);
}

inline SpinDescriptor<symm_traits::SU2Tag> operator-(SpinDescriptor<symm_traits::SU2Tag> rhs)
{
    rhs.diff_ = -rhs.diff_;
    return rhs;
}

inline std::ostream & operator<<(std::ostream & os, SpinDescriptor<symm_traits::SU2Tag> rhs)
{
    os << "Spin: " << rhs.get() << ", delta: " << rhs.action();
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
