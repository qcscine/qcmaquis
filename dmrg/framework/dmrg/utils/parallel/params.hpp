/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef PARALLEL_PARAMS_HPP
#define PARALLEL_PARAMS_HPP

namespace parallel {

    extern int groups_granularity;

    #ifdef USE_AMBIENT
    class parameters {
    public:
        static parameters& instance();
        void init();
    private:
        parameters(parameters const&);
        parameters& operator=(parameters const&);
        parameters(){}

        bool isset(const char* env);
        int getint(const char* env);
    };

    inline parameters& parameters::instance(){
        static parameters singleton;
        return singleton;
    }
    inline void parameters::init(){
        groups_granularity = !isset("DMRG_GROUPS_GRANULARITY") ? 1 
                           : getint("DMRG_GROUPS_GRANULARITY");

        std::cout << "Groups granularity: " << groups_granularity << "\n";
    }
    inline bool parameters::isset(const char* env){
        return (std::getenv( env ) != NULL);
    }
    inline int parameters::getint(const char* env){
        return std::atoi(std::getenv( env ));
    }
    #else
    class parameters {
    public:
        static parameters& instance(){
            static parameters singleton;
            return singleton;
        }
        void init(){
            groups_granularity = 1; 
        }
    private:
        parameters(parameters const&);
        parameters& operator=(parameters const&);
        parameters(){}
    };
    #endif

    extern parameters& params;
}

#endif
