#ifndef DMRGPARAMETERS_H
#define DMRGPARAMETERS_H

#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

#include <boost/any.hpp>

#include <fstream>
#include <iostream>

namespace conversion
{
    // this can be specialized to provide conversion for types that cannot be read
    // with the boost::any_cast, or whatever program options uses for the as<>
    // method
    template<class T> struct get_
    {
        T operator()(boost::program_options::variables_map & vm_, std::string const & key)
        {
            try {
                return vm_[key].as<T>();
            } catch (std::exception &e) {
                std::cerr << "Exception raised reading parameter " << key << " to type " << typeid(T).name() << std::endl;
                throw e;
            }
        }
    };
    
    template<class T> struct get_<std::vector<T> >
    {
        std::vector<T> operator()(boost::program_options::variables_map& vm_, std::string& key)
        {
            std::string raw = vm_[key].as<std::string>();
            std::vector<T> ret;
            
            typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
            boost::char_separator<char> sep(",-");
            tokenizer tokens(raw, sep);
            std::transform(tokens.begin(), tokens.end(), std::back_inserter(ret),
                           boost::lexical_cast<T, std::string>);
            return ret;
        }
    };
}

class BaseParameters
{
public:
    template<class T> T get(std::string const & key)
    {
        if (sv.count(key) > 0) { // has been changed using set()
            return boost::any_cast<T>(sv[key]);
        } else {
            conversion::get_<T> g;
            return g(vm, key);
        }
    }
    
    template<class T> void set(std::string const & key, T const & value)
    {
        sv[key] = value;
    }
    
protected:
    boost::program_options::variables_map vm;
    std::map<std::string, boost::any> sv;
};

class DmrgParameters : public BaseParameters
{
public:
    DmrgParameters(std::ifstream& param_file)
    {
        using namespace boost::program_options;
        
        try {
            options_description config("Config-file options");
            config.add_options()
            
            ("truncation_initial",value<double>(),"Initial value for the truncation error")
            ("truncation_sweep_factor",value<double>()->default_value(0.5),"After each sweep, the TE is decreased by this factor")
            
            ("max_bond_dimension",value<std::size_t>(),"")
            
            ("alpha_initial",value<double>()->default_value(1e-3),"")
            ("alpha_sweep_factor",value<double>()->default_value(0.2),"")
            
            ("eigensolver",value<std::string>()->default_value(std::string("ARPACK")),"")
            
            ("nsweeps",value<int>(),"")
            ;
            
            store(parse_config_file(param_file, config), vm);
            notify(vm);
        } catch(std::exception &e) {
            std::cerr << "Error reading parameter file." << std::endl;
            throw e;
        }
    }
};
    
class ModelParameters : public BaseParameters
{
public:
    ModelParameters(std::ifstream& param_file)
    {
        using namespace boost::program_options;
        
        try {
            options_description config("Config-file options");
            config.add_options()
            
            ("model", value<std::string>(), "")
            ("lattice", value<std::string>(), "")
            ("dof", value<std::string>(), "")
            
            ("L", value<int>(), "")
            ("W", value<int>(), "")
            
            ("theta", value<double>(), "")
            ("Jxy", value<double>(), "")
            ("Jz", value<double>(), "")
            ;
            
            store(parse_config_file(param_file, config), vm);
            notify(vm);
        } catch(std::exception &e) {
            std::cerr << "Error reading parameter file." << std::endl;
            throw e;
        }
    }
};

#endif
