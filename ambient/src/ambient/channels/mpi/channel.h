#ifndef AMBIENT_CHANNELS_MPI_CHANNEL
#define AMBIENT_CHANNELS_MPI_CHANNEL
#include "ambient/utils/singleton.hpp"
#include "ambient/channels/mpi/packets/packet_t.h"
#include "ambient/channels/mpi/packets/types.h"
#include "ambient/channels/mpi/packets/packet.h"
#include "ambient/channels/mpi/groups/group.h"
#include "ambient/channels/mpi/groups/multirank.h"
#include "ambient/utils/delegate.hpp"

namespace ambient { namespace channels { namespace mpi {

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

        group* world();
        void  init();
        void  finalize();
        std::pair<size_t*,size_t> id();
        packet_t& get_block_packet_type(size_t len);
        pipe* add_pipe(const packet_t& type, pipe::direction flow);
        pipe* get_pipe(const packet_t& type, pipe::direction flow);
        void  add_handler(const packet_t& type, void(*callback)(packet&));
        size_t get_volume() const;
        void  emit(packet* pack);
        void  spin();
        ~channel();

        group* ambient;
        std::list<pipe*> qs;
        pthread_t thread;
        pthread_mutex_t mutex;
        bool active;
    };

} } }

namespace ambient {
    extern channels::mpi::channel& channel;
}

#include "ambient/channels/mpi/channel.hpp"
#include "ambient/channels/mpi/groups/group.hpp"
#include "ambient/channels/mpi/packets/packet.hpp"
#include "ambient/channels/mpi/packets/packet_t.hpp"
#endif
