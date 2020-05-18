/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
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
#include "dmrg/utils/parallel/tracking.hpp"
#include "dmrg/utils/parallel.hpp"

// MPI stuff
#include <utils/maquis_mpi.h>

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
template<class Matrix, class SymmGroup> class SiteOperator;

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
#ifdef USE_AMBIENT
namespace ambient { inline namespace numeric {
    template <typename Matrix, int IB> class tiles;
    template <typename T, class Allocator> class matrix;
} }
namespace storage {
    template<typename T> 
    struct constrained<ambient::tiles<ambient::matrix<T> > > {
        typedef ambient::tiles<ambient::matrix<T> > type;
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
            size_t loop_max = o.aux_dim();
            for(size_t b = 0; b < loop_max; ++b){
                assert( o[b].reasonable() );
                for (std::size_t k = 0; k < o[b].n_blocks(); ++k){
                    Matrix& m = o[b][k];
#ifdef USE_AMBIENT
                    for(int j = 0; j < m.nt; ++j)
                    for(int i = 0; i < m.mt; ++i){
                        if(ambient::weak(m.tile(i,j))) continue;
                        if(ambient::ext::naked(m.tile(i,j)).state != ambient::locality::local) continue;
                        char* data = (char*)ambient::ext::naked(m.tile(i,j));
                        ofs.write(data, m.tile(i,j).num_cols() * m.tile(i,j).num_rows() *
                                  sizeof(typename Matrix::value_type)/sizeof(char));
                        std::free(data);
                        ambient::ext::naked(m.tile(i,j)).data = NULL;
                    }
#else
                    for (std::size_t c = 0; c < num_cols(m); ++c)
                        ofs.write((char*)(&m(0, c)), num_rows(m)*
                                 sizeof(typename Matrix::value_type)/sizeof(char));
                    m = Matrix();
#endif
                }
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
            size_t loop_max = o.aux_dim();
            for(size_t b = 0; b < loop_max; ++b){
                for (std::size_t k = 0; k < o[b].n_blocks(); ++k){
#ifdef USE_AMBIENT
                    Matrix& m = o[b][k];
                    for(int j = 0; j < m.nt; ++j)
                    for(int i = 0; i < m.mt; ++i){
                        if(ambient::weak(m.tile(i,j))) continue;
                        if(ambient::ext::naked(m.tile(i,j)).state != ambient::locality::local) continue;
                        ambient::ext::naked(m.tile(i,j)).data = std::malloc(ambient::ext::naked(m.tile(i,j)).spec.extent);
                        ifs.read((char*)ambient::ext::naked(m.tile(i,j)), m.tile(i,j).num_cols() * m.tile(i,j).num_rows() *
                                 sizeof(typename Matrix::value_type)/sizeof(char));
                    }
#else
                    o[b][k] = Matrix(o[b].left_basis()[k].second,
                                     o[b].right_basis()[k].second);
                    Matrix& m = o[b][k];
                    ifs.read((char*)(&m(0,0)), num_cols(m)*num_rows(m)*
                             sizeof(typename Matrix::value_type)/sizeof(char));
#endif
                }
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
            for (std::size_t b = 0; b < o.aux_dim(); ++b)
            for (std::size_t k = 0; k < o[b].n_blocks(); ++k){
#ifdef USE_AMBIENT
                    Matrix& m = o[b][k];
                    for(int j = 0; j < m.nt; ++j)
                    for(int i = 0; i < m.mt; ++i){
                        if(ambient::weak(m.tile(i,j))) continue;
                        if(ambient::ext::naked(m.tile(i,j)).state != ambient::locality::local) continue;
                        char* data = (char*)ambient::ext::naked(m.tile(i,j)); std::free(data);
                        ambient::ext::naked(m.tile(i,j)).data = NULL;
                    }
#else
                    o[b][k] = Matrix();
#endif
            }
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
           ~descriptor(){
                this->join();
            }
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
            enum { core, storing, uncore, prefetching } state;
            bool dumped;
            size_t sid;
            boost::thread* worker;
            size_t record;
        };

        template<class T> class serializable : public descriptor {
        public: 
            ~serializable(){
                if (dumped) std::remove(disk::fp(sid).c_str()); // only delete existing file, too slow otherwise on NFS or similar
            }
            serializable& operator = (const serializable& rhs){
                this->join();
                if(dumped) std::remove(disk::fp(sid).c_str());
                descriptor::operator=(rhs);
                return *this;
            }
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
                    state = storing;
                    dumped = true;
                    parallel::sync();
                    this->thread(new boost::thread(evict_request<T>(disk::fp(sid), (T*)this)));
                }
                assert(this->state != prefetching); // evict of prefetched
            }
            void drop(){
                if(dumped) std::remove(disk::fp(sid).c_str());
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
            maquis::cout << " MPI process # "<< maquis::mpi__->getGlobalRank() << ": temporary storage enabled in " << path << "\n";
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

#ifdef USE_AMBIENT
    template<class Matrix, class SymmGroup, class Scheduler = parallel::scheduler_nop> 
    static void migrate(const block_matrix<Matrix, SymmGroup>& tc, const Scheduler& scheduler = Scheduler()){
        block_matrix<Matrix, SymmGroup>& t = const_cast<block_matrix<Matrix, SymmGroup>&>(tc);
        for(int i = 0; i < t.n_blocks(); ++i){
            parallel::guard proc(scheduler(i));
            ambient::migrate(t[i]);
        }
    }
    template<class Matrix, class SymmGroup, class Scheduler = parallel::scheduler_nop> 
    static void migrate(const SiteOperator<Matrix, SymmGroup>& tc, const Scheduler& scheduler = Scheduler()){
        SiteOperator<Matrix, SymmGroup>& t = const_cast<SiteOperator<Matrix, SymmGroup>&>(tc);
        for(int i = 0; i < t.n_blocks(); ++i){
            parallel::guard proc(scheduler(i));
            ambient::migrate(t[i]);
        }
    }
    template<class Matrix, class SymmGroup, class Scheduler = parallel::scheduler_nop> 
    static void migrate(const MPSTensor<Matrix, SymmGroup>& tc, const Scheduler& scheduler = Scheduler()){
        migrate(tc.data(), scheduler); 
    }
#else
    template<typename T, class Scheduler>
    static void migrate(T& t, const Scheduler& scheduler){ }

    template<typename T>
    static void migrate(T& t)
    {
        migrate(t, parallel::scheduler_nop());
    }
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

#endif
