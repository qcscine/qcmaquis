/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef STORAGE_ARCHIVE_H
#define STORAGE_ARCHIVE_H

#include "dmrg/utils/parallel.hpp"
#include <boost/utility.hpp>
#include <alps/hdf5.hpp>
#include <alps/utility/encode.hpp>

namespace storage {

    inline std::string once(std::string fp){
        if(!parallel::uniq()) return fp+"."+parallel::rank_str();
        return fp;
    }

    inline void uniq(std::string fp){
        if(!parallel::uniq()) std::remove(once(fp).c_str());
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

        std::vector<std::string> list_children(const std::string& path) const
        {
            return impl->list_children(path);
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
