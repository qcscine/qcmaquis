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
    void clear() { }
    int get() const { return 0; }

    bool operator==(SpinDescriptor const & rhs) const { return true; }
    bool operator!=(SpinDescriptor const & rhs) const { return false; }

    template<class Archive>
    void save(Archive & ar) const
    {}

    template<class Archive>
    void load(Archive & ar)
    {}

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

    SpinDescriptor() : twoS(0), twoSaction(0) {}
    SpinDescriptor(spin_t twoS_) : twoS(twoS_), twoSaction(0) {}
    SpinDescriptor(spin_t twoS_, spin_t twoSaction_, spin_t output = 0) : twoS(twoS_), twoSaction(twoSaction_) {}

    SpinDescriptor & operator += (SpinDescriptor rhs)
    {
        // apply action of operator rhs 
        twoS += rhs.twoSaction;
        return *this;
    }

    spin_t get() const { return twoS; }
    spin_t action() const { return twoSaction; }
    spin_t input() const { return input_; }
    spin_t output() const { return output_; }

    void clear()
    {
        twoS = 0; twoSaction = 0;
    }

    bool operator==(SpinDescriptor const & rhs) const { return (twoS == rhs.twoS && twoSaction == rhs.twoSaction); }
    bool operator!=(SpinDescriptor const & rhs) const { return !(*this==rhs); }
    friend SpinDescriptor operator-(SpinDescriptor);
    friend std::ostream & operator<<(std::ostream & os, SpinDescriptor);

    template<class Archive>
    void save(Archive & ar) const
    {
        ar["twoS"] << twoS;
        ar["twoSaction"] << twoSaction;
    }

    template<class Archive>
    void load(Archive & ar)
    {
        ar["twoS"] >> twoS;
        ar["twoSaction"] >> twoSaction;
    }

    template <class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & twoS & twoSaction;
    }

private:
    spin_t twoS;
    spin_t twoSaction; // only used for operators in the MPO
    spin_t input_;
    spin_t output_;
};

// Attention: not symmetric
inline SpinDescriptor<symm_traits::SU2Tag> couple(SpinDescriptor<symm_traits::SU2Tag> a, SpinDescriptor<symm_traits::SU2Tag> b)
{
    return a += b;
}

inline SpinDescriptor<symm_traits::SU2Tag> operator-(SpinDescriptor<symm_traits::SU2Tag> rhs)
{
    rhs.twoSaction = -rhs.twoSaction;
    return rhs;
}

inline std::ostream & operator<<(std::ostream & os, SpinDescriptor<symm_traits::SU2Tag> rhs)
{
    os << "Spin: " << rhs.get() << ", Spin action: " << rhs.twoSaction;
    return os;
}

#endif
