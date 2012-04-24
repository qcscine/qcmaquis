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


namespace conversion
{
    // this can be specialized to provide conversion for types that cannot be read
    // with the boost::any_cast, or whatever program options uses for the as<>
    // method
    template<class T> struct get_
    {
        T operator()(std::string const & val)
        {
            try {
                return boost::lexical_cast<T>(val);
            } catch (std::exception &e) {
                std::cerr << "Exception raised casting " << val << " to type " << typeid(T).name() << std::endl;
                throw e;
            }
        }
    };
    
	// eliminating quatation marks around strings
	template <> struct get_<std::string>
    {
		std::string operator()(std::string const & val)
        {
            std::string ret = val;
			boost::trim_if(ret, boost::is_any_of("\"'"));
            return ret;
        }
    };	
    
    template<class T> struct get_<std::vector<T> >
    {
        std::vector<T> operator()(std::string const & val)
        {
            //            cerr << "reading " << key << std::endl;
            std::string raw = val;
			boost::trim_if(raw, boost::is_any_of("\"'"));
            std::vector<T> ret;
            
            typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
            boost::char_separator<char> sep(",");
            tokenizer tokens(raw, sep);
            std::transform(tokens.begin(), tokens.end(), std::back_inserter(ret),
                           boost::lexical_cast<T, std::string>);
            return ret;
        }
    };
}

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
            std::cerr << "Exception thrown when parsing parameters:" << std::endl;
            std::cerr << e.what() << std::endl;
            exit(1);
        }
    }
    
    bool is_set (std::string const & key)
    {
        return defined(key);
    }
    
    template<class T> T get(std::string const & key)
    {
        if (!defined(key))
            if (defaults.count(key) > 0)
                (*this)[key] = defaults[key];
            else
                boost::throw_exception(std::runtime_error("parameter " + key + " not defined"));
        conversion::get_<T> g;
        return g((*this)[key]);
    }
    
    template<class T> void set(std::string const & key, T const & value)
    {
        (*this)[key] = boost::lexical_cast<std::string>(value);
    }
    
    BaseParameters get_at_index(std::string const & var, std::size_t val, int* counter = NULL)
    {
        BaseParameters p(*this);
        
        boost::regex expression("^(.*)\\[" + var + "\\]$");
        boost::smatch what;
        for (alps::Parameters::const_iterator it=p.begin();it != p.end();++it) {
            std::string key = it->key();
            if (boost::regex_match(key, what, expression)) {
                conversion::get_<std::vector<value_type> > g;
                std::vector<value_type> v = g((*this)[key]);
                if (val < v.size())
                    p[what.str(1)] = v[val];
                else
                    p[what.str(1)] = *(v.rbegin());
                
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
