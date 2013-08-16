/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *               2011-2012 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef STORAGE_H
#define STORAGE_H

#include <boost/shared_ptr.hpp>
#include <boost/thread.hpp>
#include <boost/filesystem.hpp>

#include <iostream>
#include <fstream>

#include "utils.hpp"
#include "utils/timings.h"

#include "dmrg/utils/BaseParameters.h"
#include "dmrg/utils/tracking.h"

#ifdef HAVE_ALPS_HDF5
#include "dmrg/utils/archive.h"
#include "dmrg/utils/logger.h"
namespace storage {
    extern Logger<storage::archive> log;
}
#endif

template<class Matrix, class SymmGroup> class Boundary;
template<class Matrix, class SymmGroup> class MPSTensor;
template<class Matrix, class SymmGroup> class block_matrix;

namespace alps { namespace numeric {
    template <typename T, typename MemoryBlock> class matrix;
} }
namespace storage {
    template<class T> 
    struct constrained { 
        typedef T type; 
    };
    template<typename T> 
    struct constrained<alps::numeric::matrix<T, std::vector<T> > > {
        typedef alps::numeric::matrix<T, std::vector<T> > type;
    };
}
#ifdef AMBIENT
namespace ambient { namespace numeric {
    template <typename Matrix> class tiles;
    template <typename T, class Allocator> class matrix;
} }
namespace storage {
    template<typename T> 
    struct constrained<ambient::numeric::tiles<ambient::numeric::matrix<T, ambient::default_allocator<T> > > > {
        typedef ambient::numeric::tiles<ambient::numeric::matrix<T, ambient::constrained_allocator<T> > > type;
    };
}
#endif

namespace storage {

    class nop {
    public:
        template<class T> static void prefetch(T& o){}
        template<class T> static void fetch(T& o){}
        template<class T> static void evict(T& o){}
        template<class T> static void drop(T& o){}
        static void sync(){}
    };

#ifdef AMBIENT
    class disk : public nop {
    public:
        typedef ambient::memory::mmap impl;
        template<class T> class serializable : public impl::descriptor { };

        class pin_archive : public ambient::memory::iarchive<pin_archive> {
        public:
            template<class T, class A> void load_override(ambient::numeric::matrix<T,A>& t, int){
                ambient::safe_raw(t).spec.mmap = &pivot;
                if(ambient::safe_raw(t).state != ambient::remote)
                pivot.size += ambient::memory::aligned(ambient::safe_raw(t).spec.extent);
            }
            template<class T> void load_override(T& t, int stub){
                ambient::memory::iarchive<pin_archive>::load_override(t,stub);
            }
            pin_archive(impl::descriptor& p) : pivot(p){}
        private:
            impl::descriptor& pivot;
        };

        class load_archive : public ambient::memory::iarchive<load_archive> {
        public:
            template<class T, class A> void load_override(ambient::numeric::matrix<T,A>& t, int){
                if(ambient::safe_raw(t).data) pivot.update(ambient::safe_raw(t).data);
            }
            template<class T> void load_override(T& t, int stub){
                ambient::memory::iarchive<load_archive>::load_override(t,stub);
            }
            load_archive(impl::descriptor& p) : pivot(p){}
        private:
            impl::descriptor& pivot;
        };
       
        static void init(const std::string& path){ 
            impl::init(path);
        }
        template<class T> static void fetch(serializable<T>& t)   { if(impl::enabled() && t.remap()) load_archive(t) >> static_cast<T&>(t); }
        template<class T> static void pin(serializable<T>& t)     { if(impl::enabled()){ pin_archive(t) >> static_cast<T&>(t); t.reserve(); } }
        template<class T> static void prefetch(serializable<T>& t){ if(impl::enabled()) t.prefetch(); }
        template<class T> static void evict(serializable<T>& t)   { if(impl::enabled()) t.unmap(); }
        template<class T> static void drop(serializable<T>& t)    { if(impl::enabled()) t.drop();  }

        template<class Matrix, class SymmGroup> 
        static void evict(MPSTensor<Matrix, SymmGroup>& t){
            if(!ambient::channel.db_dim()) return;
            ambient::scope<ambient::dedicated> i;
            for(int i = 0; i < t.data().n_blocks(); ++i) 
            ambient::migrate(t.data()[i]);
        }
    };

#else

    template<class T> class evict_request {};
    template<class T> class fetch_request {};
    template<class T> class drop_request {};

    template<class Matrix, class SymmGroup>
    class evict_request< Boundary<Matrix, SymmGroup> > {
    public:
        evict_request(std::string fp, Boundary<Matrix, SymmGroup>* ptr) : fp(fp), ptr(ptr) { }
        void operator()(){
            std::ofstream ofs(fp.c_str(), std::ofstream::binary);
            Boundary<Matrix, SymmGroup>& o = *ptr;
            for (std::size_t i = 0; i < o.aux_dim(); ++i)
                for (std::size_t k = 0; k < o[i].n_blocks(); ++k){
                    for (std::size_t c = 0; c < num_cols(o[i][k]); ++c)
                        ofs.write((char*)(&o[i][k](0, c)),
                                 num_rows(o[i][k]) *
                                 sizeof(typename Matrix::value_type)/sizeof(char));
                    o[i][k] = Matrix();
                }
            ofs.close();
        }
    private:
        std::string fp;
        Boundary<Matrix, SymmGroup>* ptr;
    };

    template<class Matrix, class SymmGroup>
    class fetch_request< Boundary<Matrix, SymmGroup> > {
    public:
        fetch_request(std::string fp, Boundary<Matrix, SymmGroup>* ptr) : fp(fp), ptr(ptr) { }
        void operator()(){
            std::ifstream ifs(fp.c_str(), std::ifstream::binary);
            Boundary<Matrix, SymmGroup>& o = *ptr;
            for (std::size_t i = 0; i < o.aux_dim(); ++i)
                for (std::size_t k = 0; k < o[i].n_blocks(); ++k){
                    o[i][k] = Matrix(o[i].left_basis()[k].second,
                                     o[i].right_basis()[k].second);
                    ifs.read((char*)(&o[i][k](0,0)),
                             o[i].left_basis()[k].second*
                             o[i].right_basis()[k].second*
                             sizeof(typename Matrix::value_type)/sizeof(char));
                }
            ifs.close();
        }
    private:
        std::string fp;
        Boundary<Matrix, SymmGroup>* ptr;
    };

    template<class Matrix, class SymmGroup>
    class drop_request< Boundary<Matrix, SymmGroup> > {
    public:
        drop_request(std::string fp, Boundary<Matrix, SymmGroup>* ptr) : fp(fp), ptr(ptr) { }
        void operator()(){
            Boundary<Matrix, SymmGroup>& o = *ptr;
            for (std::size_t i = 0; i < o.aux_dim(); ++i)
            for (std::size_t k = 0; k < o[i].n_blocks(); ++k)
            o[i][k] = Matrix();
        }
    private:
        std::string fp;
        Boundary<Matrix, SymmGroup>* ptr;
    };

    class disk : public nop {
    public:
        class descriptor {
        public:
            descriptor() : state(core), dumped(false), sid(disk::index()), worker(NULL) {}
            void thread(boost::thread* t){
                this->worker = t;
                disk::track(this);
            }
            void join(){
                if(this->worker){
                    this->worker->join();
                    delete this->worker;
                    this->worker = NULL;
                    disk::untrack(this);
                }
            }
            boost::thread* worker;
            enum { core, storing, uncore, prefetching } state;
            size_t record;
            bool dumped;
            size_t sid;
        };

        template<class T> class serializable : public descriptor {
        public: 
            void fetch(){
                if(this->state == core) return;
                else if(this->state == prefetching) this->join();
                assert(this->state != storing); // isn't prefetched prior load
                assert(this->state != uncore);  // isn't prefetched prior load
                this->state = core;
            }
            void prefetch(){
                if(this->state == core) return;
                else if(this->state == prefetching) return;
                else if(this->state == storing) this->join();

                state = prefetching;
                this->thread(new boost::thread(fetch_request<T>(disk::fp(sid), (T*)this)));
            }
            void evict(){
                if(state == core){
                    if(!dumped){
                        state = storing;
                        dumped = true;
                        this->thread(new boost::thread(evict_request<T>(disk::fp(sid), (T*)this)));
                    }else{
                        state = uncore;
                        drop_request<T>(disk::fp(sid), (T*)this)();
                    }
                } 
                assert(this->state != prefetching); // evict of prefetched
            }
            void drop(){
                if(state == core) drop_request<T>(disk::fp(sid), (T*)this)();
                assert(this->state != storing);     // drop of already stored data
                assert(this->state != uncore);      // drop of already stored data
                assert(this->state != prefetching); // drop of prefetched data
            }
        };

        static disk& instance(){
            static disk singleton;
            return singleton;
        }
        static void init(const std::string& path){
            maquis::cout << "Temporary storage enabled in " << path << "\n";
            instance().active = true;
            instance().path = path;
        }
        static bool enabled(){
            return instance().active;
        }
        static std::string fp(size_t sid){
            return (instance().path + boost::lexical_cast<std::string>(sid));
        }
        static size_t index(){
            return instance().sid++;
        }
        static void track(descriptor* d){ 
            d->record = instance().queue.size();
            instance().queue.push_back(d);
        }
        static void untrack(descriptor* d){ 
            instance().queue[d->record] = NULL;
        }
        static void sync(){
            for(int i = 0; i < instance().queue.size(); ++i)
                if(instance().queue[i]) instance().queue[i]->join();
            instance().queue.clear();
        }
        template<class T> static void fetch(serializable<T>& t)   { if(enabled()) t.fetch();    }
        template<class T> static void prefetch(serializable<T>& t){ if(enabled()) t.prefetch(); }
        template<class T> static void evict(serializable<T>& t)   { if(enabled()) t.evict();    }
        template<class T> static void drop(serializable<T>& t)    { if(enabled()) t.drop();     }
        template<class T> static void pin(serializable<T>& t)     { }

        template<class Matrix, class SymmGroup> 
        static void evict(MPSTensor<Matrix, SymmGroup>& t){ }

        disk() : active(false), sid(0) {}
        std::vector<descriptor*> queue;
        std::string path;
        bool active; 
        size_t sid;
    };
#endif

    inline static void setup(BaseParameters& parms){
        if(!parms["storagedir"].empty()){
            boost::filesystem::path dp = boost::filesystem::unique_path(parms["storagedir"].as<std::string>() + std::string("/storage_temp_%%%%%%%%%%%%/"));
            try {
                boost::filesystem::create_directories(dp);
            } catch (...) {
                maquis::cerr << "Error creating dir/file at " << dp << ". Try different 'storagedir'.\n";
                throw;
            }
            storage::disk::init(dp.string());
        }else{
            maquis::cout << "Temporary storage is disabled\n";
        }
    }
}

#ifdef AMBIENT
BOOST_SERIALIZATION_REGISTER_ARCHIVE(storage::disk::pin_archive)
BOOST_SERIALIZATION_REGISTER_ARCHIVE(storage::disk::load_archive)
#endif
#endif
