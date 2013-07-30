/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef DMRG_UTILS_PARAMETER_PROXY_H
#define DMRG_UTILS_PARAMETER_PROXY_H

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

namespace parameters {
    
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
                    maquis::cerr << "Exception raised casting " << val << " to type " << typeid(T).name() << std::endl;
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
                std::string raw = val;
                boost::trim_if(raw, boost::is_any_of("\"'"));
                std::vector<T> ret;
            
                typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
                boost::char_separator<char> sep(",");
                tokenizer tokens(raw, sep);
                std::transform(tokens.begin(), tokens.end(), std::back_inserter(ret),
                               static_cast<T (*)(std::string const&)>(boost::lexical_cast<T, std::string>));
                return ret;
            }
        };
    }
    
    
    #define FOREACH_PROXY_NUMERIC_TYPE(CALLBACK)            \
    CALLBACK(double)                                        \
    CALLBACK(int)                                           \
    CALLBACK(bool)                                          \
    CALLBACK(float)                                         \
    CALLBACK(long)                                          \
    CALLBACK(std::size_t)

    #define FOREACH_PROXY_STRING_TYPE(CALLBACK)            \
    CALLBACK(std::string)

    class proxy {
    public:
        proxy(std::string & val) : val_(val) { }
        
        template <typename T>
        T as() const
        {
            conversion::get_<T> g;
            return g(val_);
        }
        
        template <typename T>
        operator T () const
        {
            return as<T>();
        }
        
        bool operator == (proxy const& rhs) const
        {
            return ( val_ == rhs.val_ );
        }
        
        bool operator == (const char* rhs) const
        {
            return ( val_ == std::string(rhs) );
        }
        
        bool operator != (const char* rhs) const
        {
            return ( val_ != std::string(rhs) );
        }
        
        template <typename T>
        proxy& operator= (T const& new_val)
        {
            val_ = boost::lexical_cast<std::string>(new_val);
            return *this;
        }
        
        std::string const& str() const
        {
            return val_;
        }
        
        const char* c_str() const
        {
            return val_.c_str();
        }
        
        bool empty() const
        {
            return val_.empty();
        }
        
        
        #define PROXY_BINARY_OP_DECL(U,OP,T)         \
        friend                                       \
        U operator OP (T lhs, proxy const& rhs)      \
        {                                            \
        return lhs OP rhs.as<T>();                   \
        }                                            \
        friend                                       \
        U operator OP (proxy const& rhs, T lhs)      \
        {                                            \
        return rhs.as<T>() OP lhs;                   \
        }
        
        #define PROXY_NUM_OPERATORS_DECL(T)            \
        PROXY_BINARY_OP_DECL(T,+,T)            \
        PROXY_BINARY_OP_DECL(T,*,T)            \
        PROXY_BINARY_OP_DECL(T,-,T)            \
        PROXY_BINARY_OP_DECL(T,/,T)            \
        PROXY_BINARY_OP_DECL(bool,>,T)         \
        PROXY_BINARY_OP_DECL(bool,>=,T)        \
        PROXY_BINARY_OP_DECL(bool,<,T)         \
        PROXY_BINARY_OP_DECL(bool,<=,T)

        #define PROXY_COMP_OPERATORS_DECL(T)       \
        PROXY_BINARY_OP_DECL(bool,==,T)            \
        PROXY_BINARY_OP_DECL(bool,!=,T)

        FOREACH_PROXY_NUMERIC_TYPE(PROXY_NUM_OPERATORS_DECL)
        FOREACH_PROXY_NUMERIC_TYPE(PROXY_COMP_OPERATORS_DECL)
        FOREACH_PROXY_STRING_TYPE(PROXY_COMP_OPERATORS_DECL)
        
        #undef PROXY_BINARY_OP_DECL
        #undef PROXY_NUM_OPERATORS_DECL
        #undef PROXY_COMP_OPERATORS_DECL
        #undef FOREACH_PROXY_NUMERIC_TYPE
        #undef FOREACH_PROXY_STRING_TYPE

        friend std::ostream& operator<< (std::ostream& os, proxy const& p)
        {
            os << p.val_;
            return os;
        }
            
    private:
        std::string & val_;
    };
}

#endif
