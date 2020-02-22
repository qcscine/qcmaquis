/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *               2019 by Leon Freitag <lefreita@ethz.ch>
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

#if !defined(BASEPARAMETERS_H) && !defined(DMRGPARAMETERS_H)
#define BASEPARAMETERS_H

#include "utils/io.hpp"

#include <string>
#include <fstream>
#include <iostream>
#include <memory>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/split_member.hpp>

#include "dmrg/utils/parameter_proxy.h"

namespace parameters {
    class value {
    public:
        value () : val_(""), empty_(true) { }

        template <class T>
        value (const T & val)
        : val_(boost::lexical_cast<std::string>(val))
        , empty_(false)
        { }

        std::string get() const {return val_;}
        bool empty() const {return empty_;}
    private:
        std::string val_;
        bool empty_;

    };
}

class BaseParameters
{
    public:

        typedef std::map<std::string, std::string> map_t;
        typedef typename map_t::value_type value_type;

        BaseParameters();
        BaseParameters(const BaseParameters & p);
        BaseParameters(std::istream& param_file);

        ~BaseParameters();

        // BaseParameters(alps::Parameters const& p);

        std::list<value_type> get_range() const;

        template <class Stream>
        void print_description(Stream& os) const;

        bool is_set (std::string const & key) const;
        bool defined (std::string const & key) const;

        parameters::proxy operator[](std::string const& key);

        template<class T> T get(std::string const & key);

        template<class T> void set(std::string const & key, T const & value);
        void set(std::string const & key, const char value[]);
        void erase(std::string const & key);

        BaseParameters iteration_params(std::string const & var, std::size_t val);

        BaseParameters & operator<<(BaseParameters const& p);

        friend void swap(BaseParameters& lhs, BaseParameters& rhs)
        {
            using std::swap;
            swap(lhs.defaults, rhs.defaults);
            swap(lhs.descriptions, rhs.descriptions);
            lhs.impl_.swap(rhs.impl_);
        }

        BaseParameters& operator=(BaseParameters rhs);
        BaseParameters& operator=(BaseParameters&& rhs);

        // for Boost::serialization
        BOOST_SERIALIZATION_SPLIT_MEMBER()
        template<class Archive>
        void load(Archive& ar);

        template<class Archive>
        void save(Archive& ar) const;


    protected:
        friend class boost::serialization::access;

        friend std::ostream& operator<<(std::ostream& os, const BaseParameters& p);

        void add_option(std::string const & name,
                        std::string const & desc,
                        parameters::value const & val = parameters::value());

        map_t defaults;
        map_t descriptions;

        class Impl;
        std::unique_ptr<Impl> impl_;
};


#endif
