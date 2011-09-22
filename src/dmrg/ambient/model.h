#ifndef AMBIENT_MODEL_H
#define AMBIENT_MODEL_H
#include "ambient/ambient.h"

namespace ambient{

    class i_model {
    public:
        virtual void add_object(unsigned int* hash, unsigned int hash_len, p_object* object) = 0;
        virtual void update_object(groups::group* scope, p_object* object) = 0;
        virtual p_object* get_object(unsigned int* hash, unsigned int hash_len, unsigned int id) const = 0;

        delegate object_modified;
        delegate object_released;
    };
    class i_channel {
    public:

    };
    class i_controller {
    public:
        virtual void acquire(i_channel* channel) = 0;

    };

    class data_model : public i_model {
    private: 
        data_model();                               // constructor is private
        data_model(data_model const&);              // copy constructor is private
        data_model& operator=(data_model const&);   // assignment operator is private
    public:
        static data_model& instance();
        virtual p_object* get_object(unsigned int* hash, unsigned int hash_len, unsigned int id) const;
        virtual void add_object(unsigned int* hash, unsigned int hash_len, p_object* object);
        virtual void update_object(groups::group* scope, p_object* object);

    private:
        hash_map map;
    };

    class d_controller : public i_controller {
    private: 
        d_controller();                                 // constructor is private
        d_controller(d_controller const&);              // copy constructor is private
        d_controller& operator=(d_controller const&);   // assignment operator is private
    public:
        static d_controller& instance();
        void acquire(i_channel* channel);
    private:
        i_model* model;
        i_channel* channel;
    };

    class mpi_channel : public i_channel {
    private: 
        mpi_channel();                                // constructor is private
        mpi_channel(mpi_channel const&);              // copy constructor is private
        mpi_channel& operator=(mpi_channel const&);   // assignment operator is private
    public:
        static mpi_channel& instance();
    private:
        i_controller* controller;
    };
}
#endif
