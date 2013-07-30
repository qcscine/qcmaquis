/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#if !defined(BASEPARAMETERS_H) && !defined(DMRGPARAMETERS_H)
#define BASEPARAMETERS_H

#include <string>
#include <fstream>
#include <iostream>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>

#include <alps/parameter.h>
#include "utils/io.hpp"

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
        bool empty_;
        std::string val_;
        
    };
}
                
class BaseParameters : public alps::Parameters
{
public:
    
    BaseParameters ()
    : alps::Parameters()
    { }
    
    BaseParameters (std::istream& param_file)
    : alps::Parameters()
    {
        try {
            alps::Parameters temp(param_file);
            *static_cast<alps::Parameters*>(this) = temp;
        } catch (std::exception & e) {
            maquis::cerr << "Exception thrown when parsing parameters:" << std::endl;
            maquis::cerr << e.what() << std::endl;
            exit(1);
        }
    }
    
    bool is_set (std::string const & key)
    {
        return defined(key);
    }

    parameters::proxy operator[](std::string const& key)
    {
        if (!defined(key)) {
            if (defaults.count(key) > 0)
                alps::Parameters::operator[](key) = defaults[key];
            else
                boost::throw_exception(std::runtime_error("parameter " + key + " not defined"));
        }
        return parameters::proxy(alps::Parameters::operator[](key));
    }
    
    template<class T> T get(std::string const & key)
    {
        parameters::proxy const& p = (*this)[key];
        return p.as<T>();
    }
    
    template<class T> void set(std::string const & key, T const & value)
    {
        alps::Parameters::operator[](key) = boost::lexical_cast<std::string>(value);
    }
    
    BaseParameters get_at_index(std::string const & var, std::size_t val, int* counter = NULL)
    {
        BaseParameters p(*this);
        
        boost::regex expression("^(.*)\\[" + var + "\\]$");
        boost::smatch what;
        for (alps::Parameters::const_iterator it=p.begin();it != p.end();++it) {
            std::string key = it->key();
            if (boost::regex_match(key, what, expression)) {
                std::vector<value_type> v = (*this)[key];
                if (val < v.size())
                    p.set(what.str(1), v[val]);
                else
                    p.set(what.str(1), *(v.rbegin()));
                
                if (counter)
                    ++(*counter);
            }
        }
        
        return p;
    }
    
protected:
    void add_option(std::string const & name,
                    std::string const & desc,
                    parameters::value const & val = parameters::value())
    {
        if (!val.empty())
            defaults[name] = val.get();
        descriptions[name] = desc;        
    }
    
    std::map<std::string, std::string> defaults;
    std::map<std::string, std::string> descriptions;
    
};

#endif
