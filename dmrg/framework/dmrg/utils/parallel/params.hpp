/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2013 by Alexandr Kosenkov <alex.kosenkov@gmail.com>
 *                         by Michele Dolfi <dolfim@phys.ethz.ch>
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
