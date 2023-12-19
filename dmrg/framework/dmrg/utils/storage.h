/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

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

#ifdef HAVE_ALPS_HDF5
#include "dmrg/utils/archive.h"
#include "dmrg/utils/logger.h"
namespace storage {
    extern Logger<storage::archive> log;
}
#endif

// Forward declaration of all the serializable objects.
template<class Matrix, class SymmGroup>
class Boundary;

template<class Matrix, class SymmGroup>
class MPSTensor;

template<class Matrix, class SymmGroup>
class block_matrix;

template<class Matrix, class SymmGroup>
class SiteOperator;

// Forward declaration of the alps matrix type
namespace alps {
namespace numeric {

template <typename T, typename MemoryBlock>
class matrix;

} // namespace numeric
} // namespace alps

namespace storage {

template<class T>
struct constrained {
    typedef T type;
};

template<typename T>
struct constrained<alps::numeric::matrix<T, std::vector<T> > > {
    typedef alps::numeric::matrix<T, std::vector<T> > type;
};

} // namespace storage

namespace storage {

/**
 * @brief Base class for a storage system modality
 * Derived class can override the methods for:
 *
 *  - prefetching (start to get something from disk)
 *  - fetching (getting something from disk)
 *  - StoreToFile (to dump objects to disk)
 *  - drop (frees the memory associated to a given object)
 *  - sync (general memory synchronization)
 *
 * Later, this class is implemented only for the disk storage case.
 */
class nop {
public:
  template<class T>
  static void prefetch(T& o) {}

  template<class T>
  static void fetch(T& o) {}

  template<class T>
  static void StoreToFile(T& o) {}

  template<class T>
  static void drop(T& o) {}

  static void sync() {}
};

// The classes here above will be specialized for the boundary class.
// Note that there is no prefetch, since this just starts the fetching process.
template<class T>
class StoreToFile_request {};

template<class T>
class fetch_request {};

template<class T>
class drop_request {};

/** @brief Functor class for the StoreToFile method, specialized for the boundaries */
template<class Matrix, class SymmGroup>
class StoreToFile_request< Boundary<Matrix, SymmGroup> > {
public:
  /** @brief Class constructor */
  StoreToFile_request(std::string fp, Boundary<Matrix, SymmGroup>* ptr) : fp(fp), ptr(ptr) { }

  /** @brief Round brackets operator */
  void operator()() {
    std::ofstream ofs(fp.c_str(), std::ofstream::binary);
    Boundary<Matrix, SymmGroup>& o = *ptr;
    auto loop_max = o.aux_dim();
    for (int b = 0; b < loop_max; ++b) {
      assert( o[b].reasonable() );
      for (int k = 0; k < o[b].n_blocks(); ++k) {
        Matrix& m = o[b][k];
        for (int c = 0; c < num_cols(m); ++c)
          ofs.write((char*)(&m(0, c)), num_rows(m)*sizeof(typename Matrix::value_type)/sizeof(char));
        m = Matrix();
      }
    }
    ofs.close();
  }
private:
  std::string fp;
  Boundary<Matrix, SymmGroup>* ptr;
};

/** @brief Functor class for the fetch method, specialized for the boundaries */
template<class Matrix, class SymmGroup>
class fetch_request< Boundary<Matrix, SymmGroup> > {
public:
  /** @brief Class constructor */
  fetch_request(std::string fp, Boundary<Matrix, SymmGroup>* ptr) : fp(fp), ptr(ptr) { }

  /** @brief Round braket operator */
  void operator()() {
    std::ifstream ifs(fp.c_str(), std::ifstream::binary);
    Boundary<Matrix, SymmGroup>& o = *ptr;
    auto loop_max = o.aux_dim();
    for (int b = 0; b < loop_max; ++b) {
      for (int k = 0; k < o[b].n_blocks(); ++k) {
        o[b][k] = Matrix(o[b].left_basis()[k].second,
                         o[b].right_basis()[k].second);
        Matrix& m = o[b][k];
        ifs.read((char*)(&m(0,0)), num_cols(m)*num_rows(m)*
                 sizeof(typename Matrix::value_type)/sizeof(char));
      }
    }
    ifs.close();
  }
private:
  std::string fp;
  Boundary<Matrix, SymmGroup>* ptr;
};

/** @brief Functor class for the drop method, specialized for the boundaries */
template<class Matrix, class SymmGroup>
class drop_request< Boundary<Matrix, SymmGroup> > {
public:
  /** @brief Class constructor */
  drop_request(std::string fp, Boundary<Matrix, SymmGroup>* ptr) : fp(fp), ptr(ptr) { }

  /** @brief Round braket operator */
  void operator()() {
    Boundary<Matrix, SymmGroup>& o = *ptr;
    for (int b = 0; b < o.aux_dim(); ++b)
      for (int k = 0; k < o[b].n_blocks(); ++k)
        o[b][k] = Matrix();
  }
private:
    std::string fp;
    Boundary<Matrix, SymmGroup>* ptr;
};

/**
 * @brief Specialization of the nop class for disk storage
 * Note that this class wraps two classes:
 *
 *  - descriptor, which represents an abstract object which is associated to
 *    a given status (describing whether it's currently being fetched, stored etc)
 *    and a given boost::thread that manages it.
 *
 */
class disk : public nop {
public:

  /** @brief Descriptor class, described above */
  class descriptor {
  public:
    /** @brief Class constructor */
    descriptor() : state(core), dumped(false), sid(disk::index()), worker(NULL) {}

    /**
     * @brief Class destructor.
     * Note that the descructor calls the join method on the boost::thread -- i.e.,
     * completes the operation associated with the boost::thread.
     */
    ~descriptor() { this->join(); }

    /**
     * @brief Associated the object to a given thread.
     * Note that this method calls the disk::track method, which ensures that the corresponding
     * thread is being "followed" by the memory manager.
     */
    void thread(boost::thread* t){
      this->worker = t;
      disk::track(this);
    }

    /**
     * @brief Executes the thread.
     * Similarly to the [thread] method, also here the disk::untrack method is called under the
     * hood in order to not follow anymore this thread.
     */
    void join(){
      if(this->worker){
        this->worker->join();
        delete this->worker;
        this->worker = NULL;
        disk::untrack(this);
      }
    }

    // Enum class identifying the possible states of the object.
    enum { core, storing, uncore, prefetching } state;

    /** Class members */
    bool dumped;					    // Bool keeping track of whether the object has been written to disk.
    size_t sid;						    // Identifier of the memory
    boost::thread* worker;	  // Boost thread managing the obect
    size_t record;						// Size associated with the object.
  };

  /**
   * @brief Class representing a serializable object.
   *
   * The class is inherited by [descriptor] such that the object is "equipped" with
   * the boost::thread and with the flags describing its storing status.
   *
   * @tparam T type associated with the serialized object.
   */
  template<class T>
  class serializable : public descriptor {
  public:
    /**
     * @brief Class descructor
     *
     * Note that sid is a member of the mother class and is the unique identified associated
     * with the object to be serialized
     */
    ~serializable() {
      // only delete existing file, too slow otherwise on NFS or similar
      if (dumped)
        std::remove(disk::fp(sid).c_str());
    }

    /** @brief Copy constructor */
    serializable& operator = (const serializable& rhs){
      this->join();
      if(dumped)
        std::remove(disk::fp(sid).c_str());
      descriptor::operator=(rhs);
      return *this;
    }

    /** @brief Fetching method */
    void fetch() {
      if(this->state == core)
        return;
      else if(this->state == prefetching)
        this->join();
      // Cannot fetch something that is being stored
      assert(this->state != storing);
      assert(this->state != uncore);
      this->state = core;
    }

    /** @brief Starts the fetching */
    void prefetch() {
      // If already available on disk, does nothing. Otherwise, if it's being
      // stored, finalized the storing such that afterwards one can call "fetch".
      if(this->state == core)
        return;
      else if(this->state == prefetching)
        return;
      else if(this->state == storing)
        this->join();
      state = prefetching;
      // This thread will be joined by [fetch]. Note that here we call the
      // functor class defined above.
      this->thread(new boost::thread(fetch_request<T>(disk::fp(sid), (T*)this)));
    }

    /** @brief Storing to file method */
    void StoreToFile(){
      if(state == core) {
        state = storing;
        dumped = true;
        parallel::sync();
        this->thread(new boost::thread(StoreToFile_request<T>(disk::fp(sid), (T*)this)));
      }
      assert(this->state != prefetching);
    }

    /** @brief Free the memory */
    void drop(){
      if(dumped)
        std::remove(disk::fp(sid).c_str());
      if(state == core)
        drop_request<T>(disk::fp(sid), (T*)this)();
      assert(this->state != storing);     // drop of already stored data
      assert(this->state != uncore);      // drop of already stored data
      assert(this->state != prefetching); // drop of prefetched data
    }
  };

  /**
   * @brief Implementation of the singleton method.
   * The actual object is created here, so it's ensured to be created only once.
   */
  static disk& instance(){
    static disk singleton;
    return singleton;
  }

  /**
   * @brief Initialization of the memory.
   * @param path file where the object will be dumped.
   */
  static void init(const std::string& path){
    maquis::cout << "Temporary storage enabled in " << path << "\n";
    instance().active = true;
    instance().path = path;
  }

  /** @brief Disables the storage to file */
  static void disable() {
    instance().active = false;
    instance().path = "";
  }

  /** @brief Checks if memory dumping has been enabled */
  static bool enabled() { return instance().active; }

  /** @brief Generates the path where a given id is stored */
  static std::string fp(size_t sid){
    return (instance().path + boost::lexical_cast<std::string>(sid));
  }

  /** @brief Generates a new index, to be associated with a new object */
  static size_t index(){
    return instance().sid++;
  }

  /** @brief Adds a descriptor to the object to be tracked */
  static void track(descriptor* d){
    d->record = instance().queue.size();
    instance().queue.push_back(d);
  }

  /** @brief Removes a descriptor from the list of objects to be tracked */
  static void untrack(descriptor* d){
      instance().queue[d->record] = NULL;
  }

  /** @brief Syncs all the processes that are queued */
  static void sync(){
    for(int i = 0; i < instance().queue.size(); ++i)
      if(instance().queue[i])
        instance().queue[i]->join();
    instance().queue.clear();
  }

  /** @brief Gets a serializable object from disk */
  template<class T>
  static void fetch(serializable<T>& t) {
    if (enabled())
      t.fetch();
  }

  /** @brief Starts getting a serializable object from disk */
  template<class T>
  static void prefetch(serializable<T>& t) {
    if (enabled())
      t.prefetch();
  }

  /** @brief Stores a serializable object on disk */
  template<class T>
  static void StoreToFile(serializable<T>& t) {
    if (enabled()) {
      std::cout << "I AM ACTUALLY STORING TO FILE" << std::endl;
      t.StoreToFile();
    }
  }

  /** @brief Drops the memory associated with a serializable object */
  template<class T>
  static void drop(serializable<T>& t) {
    if (enabled())
      t.drop();
  }

  /** @brief Overload for an MPS -- never stored */
  template<class Matrix, class SymmGroup>
  static void StoreToFile(MPSTensor<Matrix, SymmGroup>& t){ }

  /** @brief Constructor for a disk (but should be made private in agreement with the singleton)*/
  disk() : active(false), sid(0) {}

  // Class members
  std::vector<descriptor*> queue; // List of objects to be tracked
  std::string path;							  // Location where the objects are stored
  bool active;									  // Whether in-disk storage is activated or not
  size_t sid;										  // Last used index
};

template<typename T, class Scheduler>
static void migrate(T& t, const Scheduler& scheduler){ }

template<typename T>
static void migrate(T& t) { migrate(t, parallel::scheduler_nop()); }

/**
 * @brief Initialization method.
 * This method chooses -- depending on the input parameters -- whether to activate in-disk storage or not.
 * Note that this method has to be called EVERY TIME a new calculation is started, otherwise the singleton
 * is not reset.
 * For this reason, setup is called in the constructor of the simulation object.
 */
inline static void setup(BaseParameters& parms)
{
  if(!parms["storagedir"].empty()) {
    auto dp = boost::filesystem::unique_path(parms["storagedir"].as<std::string>() + std::string("/storage_temp_%%%%%%%%%%%%/"));
    try {
      boost::filesystem::create_directories(dp);
    }
    catch (...) {
      maquis::cerr << "Error creating dir/file at " << dp << ". Try different 'storagedir'.\n";
      throw;
    }
    storage::disk::init(dp.string());
  }
  else {
    storage::disk::disable();
    maquis::cout << "Disk storage of boundaries disabled" << std::endl;;
  }
}

} // namespace storage

#endif
