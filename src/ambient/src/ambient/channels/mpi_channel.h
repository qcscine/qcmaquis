#ifndef AMBIENT_CHANNELS_MPI_CHANNEL
#define AMBIENT_CHANNELS_MPI_CHANNEL
#include "ambient/ambient.h"
#include "ambient/channels/ichannel.h"
#include "ambient/channels/packets/types.h"
#include "ambient/channels/packets/packet.h"
#include "ambient/channels/groups/group.h"
#include "ambient/utils/singleton.hpp"
#include "ambient/utils/delegate.hpp"
#include <list>

namespace ambient { namespace channels {

    class mpi_channel : public ichannel, public singleton< mpi_channel > 
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
            pipe(const channels::packet_t& type, direction flow);
           ~pipe();
            size_t get_bound() const;
            request* create_request();
            request* attach_request(void* memory);
            void renew(request* r);
            void send(request* r);
            void recv(request* r);
            void spin();
            const channels::packet_t& type;
            delegate packet_delivered;
            direction flow;
        private:
            std::list<request*> reqs;
            pthread_mutex_t reqs_mutex;
        };
    public:
        mpi_channel();
        static void* stream(void* instance);
        void ifetch(group* placement, size_t gid, size_t sid, size_t i, size_t j); // exclude owner?

        group* world();
        void  init();
        void  finalize();
        std::pair<size_t*,size_t> id();
        channels::packet_t& get_block_packet_type(size_t len);
        pipe* add_pipe(const channels::packet_t& type, pipe::direction flow);
        pipe* get_pipe(const channels::packet_t& type, pipe::direction flow);
        void  add_handler(const channels::packet_t& type, void(*callback)(channels::ichannel::packet&));
        size_t get_volume() const;
        void  emit(channels::ichannel::packet* pack);
        void  spin();
       ~mpi_channel();

        channels::group* ambient;
        std::list<pipe*> qs;
        pthread_t thread;
        pthread_mutex_t mutex;
        bool active;
    };

} }

#endif
