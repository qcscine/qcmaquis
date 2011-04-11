/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef DMRGPARAMETERS_H
#define DMRGPARAMETERS_H

#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include <boost/any.hpp>

#include <fstream>
#include <iostream>

#ifdef HAVE_ALPS_HDF5
#include <alps/hdf5.hpp>
#endif

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

	// eliminating quatation marks around strings
	template <> struct get_<std::string>
    {
		std::string operator()(boost::program_options::variables_map& vm_, std::string const & key)
        {
            std::string ret = vm_[key].as<std::string>();
			boost::trim_if(ret, boost::is_any_of("\"'"));
            return ret;
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

namespace detail
{
    struct SerializationHelperBase {
        virtual void operator()(alps::hdf5::oarchive & ar,
                                std::string const & path,
                                boost::program_options::variable_value const & a) const = 0;
    };
    
    template<class T>
    struct SerializationHelper : public SerializationHelperBase
    {
        void operator()(alps::hdf5::oarchive & ar,
                        std::string const & path,
                        boost::program_options::variable_value const & a) const
        {
            if (!a.empty())
                ar << alps::make_pvp(path, a.as<T>());
        }
    };
}
            

class BaseParameters
{
public:
    bool is_set (std::string const & key)
    {
        return (vm.count(key)) || (sv.count(key) > 0);
    }
    
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
    
#ifdef HAVE_ALPS_HDF5
    void serialize(alps::hdf5::oarchive & ar) const
    {
        using boost::program_options::option_description;
        using boost::shared_ptr;
        
        std::vector<shared_ptr<option_description> > opts = config.options();
        
        for (std::vector<shared_ptr<option_description> >::iterator it = opts.begin();
             it != opts.end(); ++it)
        {
            std::string name = (*it)->long_name();
            detail::SerializationHelperBase * helper = shelpers.find(name)->second.get();
            try {
                (*helper)(ar, name, vm[name]);
            } catch (std::exception &e) {
                std::cerr << "Error writing parameter " << name << std::endl;
                throw e;
            }
        }
    }
#endif
    
protected:
    boost::program_options::options_description config;
    boost::program_options::variables_map vm;
    std::map<std::string, boost::any> sv;
    std::map<std::string, boost::shared_ptr<detail::SerializationHelperBase> > shelpers;
    
    template<class T>
    void add_option(const char * name,
                    boost::program_options::typed_value<T> * val,
                    const char * desc)
    {
        config.add_options()(name, val, desc);
        shelpers[name] = boost::shared_ptr<detail::SerializationHelperBase>(new detail::SerializationHelper<T>());
    }
};

class DmrgParameters : public BaseParameters
{
public:
    DmrgParameters(std::ifstream& param_file)
    {
        using namespace boost::program_options;
        
        try {
            add_option("truncation_initial",value<double>(),"Initial value for the truncation error");
            add_option("truncation_final",value<double>(),"Final value for the truncation");
            
            add_option("init_bond_dimension",value<std::size_t>()->default_value(5),"");
            add_option("max_bond_dimension",value<std::size_t>(),"");
            
            add_option("alpha_initial",value<double>()->default_value(1e-2),"");
            add_option("alpha_final",value<double>()->default_value(1e-6),"");
            
            add_option("eigensolver",value<std::string>()->default_value(std::string("ARPACK")),"");
            add_option("arpack_tol",value<double>()->default_value(1e-8),"");
            add_option("arpack_ncv",value<int>()->default_value(20),"");
            add_option("ietl_jcd_tol",value<double>()->default_value(1e-8),"");
            add_option("ietl_jcd_gmres",value<int>()->default_value(5),"");
            add_option("ietl_jcd_maxiter",value<int>()->default_value(100),"");
            
            add_option("nsweeps",value<int>(),"");
            add_option("ngrowsweeps",value<int>(),"");
            
            add_option("resultfile",value<std::string>(),"");
            add_option("chkpfile",value<std::string>(),"");
            add_option("initfile",value<std::string>()->default_value(std::string()),"");
            
            add_option("donotsave",value<int>()->default_value(0),"");
            add_option("run_seconds",value<int>()->default_value(0),"");
            add_option("storagedir",value<std::string>()->default_value(std::string()),"");
            add_option("use_compressed",value<int>()->default_value(0),"");
            add_option("calc_h2",value<int>()->default_value(0),"");
            add_option("seed",value<int>()->default_value(42),"");
            
            add_option("init_state", value<std::string>()->default_value("default"),"");
            
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
            add_option("model", value<std::string>(), "");
            add_option("lattice", value<std::string>(), "");
            add_option("alps_lattice",value<std::string>(),"");
            
            add_option("L", value<int>(), "");
            add_option("W", value<int>(), "");
            
            add_option("Jxy", value<double>(), "");
            add_option("Jz", value<double>(), "");
            
            add_option("U", value<double>(), "");
            add_option("t", value<double>()->default_value(1.), "");
            add_option("t1", value<double>()->default_value(1.), "");
            add_option("t2", value<double>()->default_value(1.), "");
            
            add_option("theta", value<double>(), "");
            add_option("h0", value<double>()->default_value(0), "");
            add_option("pin", value<int>()->default_value(3), "");
            
            add_option("K0", value<double>()->default_value(1), "");
            add_option("K1", value<double>()->default_value(0), "");
            
            add_option("penalty", value<double>()->default_value(10), "");
            add_option("twist", value<double>()->default_value(0), "");
            add_option("move", value<double>()->default_value(0), "");
            
            add_option("u1_total_charge", value<int>()->default_value(0), "");
            add_option("u1_total_charge1", value<int>()->default_value(0), "");
            add_option("u1_total_charge2", value<int>()->default_value(0), "");
            
            store(parse_config_file(param_file, config), vm);
            notify(vm);
        } catch(std::exception &e) {
            std::cerr << "Error reading parameter file." << std::endl;
            std::cerr << e.what() << endl;
            throw e;
        }
    }
};

#endif
