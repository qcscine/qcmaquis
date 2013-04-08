#ifndef AMBIENT_CHANNELS_MPI_CHANNEL
#define AMBIENT_CHANNELS_MPI_CHANNEL
#include "ambient/utils/singleton.hpp"
#include "ambient/channels/mpi/packets/packet_t.h"
#include "ambient/channels/mpi/packets/types.h"
#include "ambient/channels/mpi/packets/packet.h"
#include "ambient/channels/mpi/groups/group.h"
#include "ambient/channels/mpi/groups/multirank.h"
#include "ambient/utils/delegate.hpp"
#include <atomic>

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

    class channel : public singleton< channel > 
    {
        class pipe {
            class request {
            public:
                request(void* memory);
                MPI_Request mpi_request;
                void* memory;
            };
        public:
            enum direction { IN, OUT, LO };
            pipe(const packet_t& type, direction flow);
            ~pipe();
            size_t get_bound() const;
            request* create_request();
            request* attach_request(void* memory);
            void renew(request* r);
            void send(request* r);
            void recv(request* r);
            void spin();
            const packet_t& type;
            delegate packet_delivered;
            direction flow;
        private:
            std::list<request*> reqs;
            pthread_mutex_t reqs_mutex;
        };
    public:
        channel();
        static void* stream(void* instance);
        void ifetch(group* placement, size_t sid, size_t x, size_t y);

        void  init();
        pipe* add_pipe(const packet_t& type, pipe::direction flow);
        pipe* get_pipe(const packet_t& type, pipe::direction flow);
        void  add_handler(const packet_t& type, void(*callback)(packet&));

        request* get(revision* r);
        request* set(revision* r, int rank);
        request* get(transformable* v);
        request* set(transformable* v, int rank);
        bool  test(request* r);

        void  emit(packet* p);
        void  spin();
       ~channel();

        group* world;
        size_t volume;
        std::list<pipe*> qs;
        std::atomic<bool> active;

        pthread_t thread;
        pthread_mutex_t mutex;
    };

} } }

namespace ambient {
    extern channels::mpi::channel& channel;
}

#include "ambient/channels/mpi/channel.hpp"
#include "ambient/channels/mpi/groups/multirank.hpp"
#include "ambient/channels/mpi/groups/group.hpp"
#include "ambient/channels/mpi/packets/packet.hpp"
#include "ambient/channels/mpi/packets/packet_t.hpp"
#endif
