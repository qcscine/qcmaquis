/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
*               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
*               2011-2013    Michele Dolfi <dolfim@phys.ethz.ch>
*               2014-2014    Sebastian Keller <sebkelle@phys.ethz.ch>
*               2019         Leon Freitag <lefreita@ethz.ch>
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
#ifndef INTEGRAL_INTERFACE_H
#define INTEGRAL_INTERFACE_H

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/complex.hpp>


namespace chem {

    // // Leon: TODO:
    // // distinguish between electronic or vibrational
    // // This way we can merge parse_integrals_vib.h here

    // // The only thing that changes in electronic vs vibrational integrals is the
    // // number of indices in the FCIDUMP file: 4 in the electronic and 6 in the vibrational
    // const int Electronic = 4, Vibrational = 6;
    //
    // template <int N = Electronic>
    // using index_type = std::array<int, N>;

    // Classes to represent integrals

    typedef std::array<int, 4> index_type;

    template <class V>
    using integral_tuple = std::pair<index_type, V>;

    template <class V>
    using integrals = std::vector<integral_tuple<V> >; // TODO: use a map later

    // Class that regulates passing of the integrals between the model and outside world.
    // The only way to feed integrals to a model is via the parameters, which ultimately
    // store strings, however if we want to interface it to other applications, we must
    // provide integrals as a reasonable data structure. This class handles the conversion
    // of an integral data structure into the string that may be fed into model parameters.
    template <class V>
    class integral_proxy
    {
        public:

            integral_proxy() = default;
            integral_proxy(const integrals<V> & ints) : ints_(ints) {};
            integral_proxy(integrals<V> && ints) : ints_(std::move(ints)) {};

            const integrals<V> & ints() { return ints_; };

        private:
            friend class boost::serialization::access;

            // Currently, we store integral indices and values in two separate vectors.
            // Indexing is done as in the FCIDUMP format, that is:
            // Orbital indices start from 1 and 2e integrals use all four indices
            // 1e integrals use the first two indices and 0,0 as the 3rd and the 4th index
            // Nuclear repulsion energy uses an index 0,0,0,0
            integrals<V> ints_;

            template <typename Archive>
            friend void serialize(Archive& ar, integral_proxy &i, const unsigned int version)
            {
                ar & i.ints_;
            }
    };


    // Serialize the structure into a string

    template <class V>
    std::string serialize(const integrals<V>& ints)
    {
        integral_proxy<V> p(ints);
        std::stringstream ss;
        boost::archive::text_oarchive oa{ss};
        oa << p;

        return ss.str();
    }
}

#endif