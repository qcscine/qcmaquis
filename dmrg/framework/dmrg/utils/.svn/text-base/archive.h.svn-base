/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *               2011-2012 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef STORAGE_ARCHIVE_H
#define STORAGE_ARCHIVE_H

#ifdef USE_AMBIENT
#include "ambient/ambient.hpp"
#endif

#include <boost/utility.hpp>
#include <alps/hdf5.hpp>
#include <alps/utility/encode.hpp>

namespace storage {

    inline std::string once(std::string fp){
        #ifdef USE_AMBIENT
        if(!ambient::master() && !ambient::scope::nested()) return fp+"."+std::to_string(ambient::rank());
        #endif
        return fp;
    }

    inline void uniq(std::string fp){
        #ifdef USE_AMBIENT
        if(!ambient::master() && !ambient::scope::nested()) std::remove(once(fp).c_str());
        #endif
    }

    class archive : boost::noncopyable {
    public:
        archive(std::string fp) : write(false), fp(fp) {
            impl = new alps::hdf5::archive(fp);
        }
        archive(std::string fp, const char* rights) : write(strcmp(rights,"w") == 0), fp(fp) {
            impl = new alps::hdf5::archive(once(fp), rights); 
        }
       ~archive(){
           delete impl;
           if(write) uniq(fp); 
        }
        bool is_group(const char* path){
            return impl->is_group(path);
        }
        bool is_scalar(const char* path){
            return impl->is_scalar(path);
        }
        bool is_data(const char* path){
            return impl->is_data(path);
        }
        template<typename T>
        void operator << (const T& obj){
            (*impl) << obj;
        }
        template<typename T>
        void operator >> (T& obj){
            (*impl) >> obj;
        }
        alps::hdf5::detail::archive_proxy<alps::hdf5::archive> operator[](std::string path){
            return (*impl)[path];
        }
    private:
        bool write;
        std::string fp;
        alps::hdf5::archive* impl;
    };
    
    inline std::string encode(std::string const & s){
        return alps::hdf5_name_encode(s);
    }
}

#endif
