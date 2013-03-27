/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Tim Ewart <timothee.ewart@gmail.com>
 *
 *****************************************************************************/

#ifndef STREAM_STORAGE_AMBIENT_H
#define STREAM_STORAGE_AMBIENT_H

#include "dmrg/block_matrix/detail/ambient.hpp"

#ifdef USE_COMPLEX
    #define dmrg_value_type std::complex<double>
#else
    #define dmrg_value_type double
#endif

template<class SymGroup>
class StreamWriteRequest_impl< Boundary<typename ambient::numeric::tiles<ambient::numeric::matrix<dmrg_value_type> >, SymGroup> >
: public StreamRequest
{
    typedef Boundary<typename ambient::numeric::tiles<ambient::numeric::matrix<dmrg_value_type> >, SymGroup> Object;
    typedef std::list< boost::tuple<void*, std::size_t, std::size_t> > ObjectList;
    
public:
    StreamWriteRequest_impl(StreamStorage * store_,
                                     Object * ptr_,
                                     boost::shared_ptr<boost::mutex> wait_mutex_)
    : store(store_)
    , wait_mutex(wait_mutex_)
    , wait_lock(*wait_mutex_)
    {
        std::size_t tag(0);
        for (std::size_t i = 0; i < ptr_->aux_dim(); ++i)
            for (std::size_t k = 0; k < (*ptr_)[i].n_blocks(); ++k) {
                typename ambient::numeric::tiles<ambient::numeric::matrix<dmrg_value_type> > tmp(num_rows((*ptr_)[i][k]), num_cols((*ptr_)[i][k]));
                ambient::numeric::save((*ptr_)[i][k],tag);
                tag += (*ptr_)[i][k].data.size();
                (*ptr_)[i][k].swap(tmp);
            }

        ambient::sync(); // the list is full
        ambient::fout.get_list(list); // get the list from io manager
    }
    
    void operator()()
    {
        std::string fp = store->master()->get_base_path() + store->object_path;
        std::ofstream of(fp.c_str(), std::ofstream::binary);

        typename ObjectList::const_iterator it = list.begin();
        for(; it != list.end(); ++it){
           of.write((char*)(dmrg_value_type*)(boost::get<0>((*it))),boost::get<2>((*it))/sizeof(char));
           free(boost::get<0>((*it))); //clean memory
        }

        of.close();
    }
    
private:
    StreamStorage * store;
    ObjectList  list;
    boost::shared_ptr<boost::mutex> wait_mutex;
    boost::lock_guard<boost::mutex> wait_lock;
};

template<class SymGroup>
class StreamReadRequest_impl< Boundary<typename ambient::numeric::tiles<ambient::numeric::matrix<dmrg_value_type> >, SymGroup> >
: public StreamRequest
{
    typedef Boundary<typename ambient::numeric::tiles<ambient::numeric::matrix<dmrg_value_type> >, SymGroup> Object;
    typedef std::list< boost::tuple<void*, std::size_t, std::size_t> > ObjectList;
    
public:
     StreamReadRequest_impl(StreamStorage * store_,
                                     Object * ptr_,
                                     boost::shared_ptr<boost::mutex> wait_mutex_)
    : store(store_)
    , wait_mutex(wait_mutex_)
    , wait_lock(*wait_mutex_)
    {
        std::size_t tag(0);
        for (std::size_t i = 0; i < ptr_->aux_dim(); ++i)
            for (std::size_t k = 0; k < (*ptr_)[i].n_blocks(); ++k) {
                ambient::numeric::load((*ptr_)[i][k],tag);
                tag += (*ptr_)[i][k].data.size();
            }

        ambient::sync(); // the list is full
        ambient::fout.get_list(list); // get the list from io manager
    }
    
    void operator()()
    {
        std::string fp = store->master()->get_base_path() + store->object_path;
        std::ifstream ifs(fp.c_str(), std::ofstream::binary);

        typename ObjectList::const_iterator it = list.begin();
        for(; it != list.end(); ++it)
           ifs.read((char*)((dmrg_value_type*)boost::get<0>((*it))),  boost::get<2>((*it))/sizeof(char));
       
        ifs.close();
    }
    
private:
    StreamStorage * store;
    ObjectList  list;
    boost::shared_ptr<boost::mutex> wait_mutex;
    boost::lock_guard<boost::mutex> wait_lock;
};

template<class SymGroup>
class StreamWriteRequest_impl< block_matrix<typename ambient::numeric::tiles<ambient::numeric::matrix<dmrg_value_type> >, SymGroup> >
: public StreamRequest
{
    typedef block_matrix<typename ambient::numeric::tiles<ambient::numeric::matrix<dmrg_value_type> >, SymGroup> Object;
    typedef std::list< boost::tuple<void*, std::size_t, std::size_t> > ObjectList;
    
public:
    StreamWriteRequest_impl(StreamStorage * store_,
                            Object * ptr_,
                            boost::shared_ptr<boost::mutex> wait_mutex_)
    : store(store_)
    , wait_mutex(wait_mutex_)
    , wait_lock(*wait_mutex_)
    {
        assert(false); // not used
        std::size_t tag(0);
        for (std::size_t k = 0; k < ptr_.n_blocks(); ++k) {
            ambient::numeric::save((*ptr_)[k],tag);
            tag +=  (*ptr_)[k].data.size();
            typename ambient::numeric::tiles<ambient::numeric::matrix<dmrg_value_type> > tmp((*ptr_)[k].num_rows(), (*ptr_)[k].num_cols());                 
            (*ptr_)[k].swap(tmp);
        }

        ambient::sync(); // the list is full
        ambient::fout.get_list(list); // get the list from io manager
    }
    
    void operator()()
    {
        std::string fp = store->master()->get_base_path() + store->object_path;
        std::ofstream of(fp.c_str(), std::ofstream::binary);

        typename ObjectList::const_iterator it = list.begin();
        for(; it != list.end(); ++it){
           of.write((char*)(dmrg_value_type*)(boost::get<0>((*it))),boost::get<2>((*it))/sizeof(char));
           free(boost::get<0>((*it))); //clean memory
        }

        of.close();
    }
    
private:
    StreamStorage * store;
    ObjectList  list;
    boost::shared_ptr<boost::mutex> wait_mutex;
    boost::lock_guard<boost::mutex> wait_lock;
};

template<class SymGroup>
class StreamReadRequest_impl< block_matrix<typename ambient::numeric::tiles<ambient::numeric::matrix<dmrg_value_type> >, SymGroup> >
: public StreamRequest
{
    typedef block_matrix<typename ambient::numeric::tiles<ambient::numeric::matrix<dmrg_value_type> >, SymGroup> Object;
    typedef std::list< boost::tuple<void*, std::size_t, std::size_t> > ObjectList;
    
public:
    StreamReadRequest_impl(StreamStorage * store_,
                            Object * ptr_,
                            boost::shared_ptr<boost::mutex> wait_mutex_)
    : store(store_)
    , wait_mutex(wait_mutex_)
    , wait_lock(*wait_mutex_)
    {
        assert(false); // not used
        std::size_t tag(0);
        for (std::size_t k = 0; k < ptr_.n_blocks(); ++k) {
            ambient::numeric::load((*ptr_)[k],tag);
            tag +=  (*ptr_)[k].data.size();
        }

        ambient::sync(); // the list is full
        ambient::fout.get_list(list); // get the list from io manager
    }
    
    void operator()()
    {
        std::string fp = store->master()->get_base_path() + store->object_path;
        std::ifstream ifs(fp.c_str(), std::ofstream::binary);

        typename ObjectList::const_iterator it = list.begin();
        for(; it != list.end(); ++it)
           ifs.read((char*)(dmrg_value_type*)(boost::get<0>((*it))),boost::get<2>((*it))/sizeof(char));

        ifs.close();
    }
    
private:
    StreamStorage * store;
    ObjectList list;
    boost::shared_ptr<boost::mutex> wait_mutex;
    boost::lock_guard<boost::mutex> wait_lock;
};

#undef dmrg_value_type

#endif
