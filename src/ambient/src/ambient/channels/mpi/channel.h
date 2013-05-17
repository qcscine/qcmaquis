#ifndef AMBIENT_CHANNELS_MPI_CHANNEL
#define AMBIENT_CHANNELS_MPI_CHANNEL
#include "ambient/utils/singleton.hpp"
#include "ambient/channels/mpi/groups/group.h"
#include "ambient/channels/mpi/groups/multirank.h"

namespace ambient { namespace channels { namespace mpi {

    using ambient::models::velvet::revision;
    using ambient::models::velvet::transformable;

    class request {
    public:
        void* operator new (size_t size){ return ambient::bulk.malloc<sizeof(request)>(); }
        void operator delete (void* ptr){ }
        request(void* memory);
        MPI_Request mpi_request;
        void* memory;
    };

    class channel : public singleton< channel > {
    public:
       ~channel();
        void  init();
        request* get(revision* r);
        request* set(revision* r, int rank);
        request* get(transformable* v);
        request* set(transformable* v, int rank);
        bool test(request* r);
    public:
        group* world;
        size_t volume;
    };

} } }

namespace ambient {
    extern channels::mpi::channel& channel;
}

#include "ambient/channels/mpi/channel.hpp"
#include "ambient/channels/mpi/groups/multirank.hpp"
#include "ambient/channels/mpi/groups/group.hpp"
#endif
