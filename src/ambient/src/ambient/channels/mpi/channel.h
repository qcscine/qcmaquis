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
                inline request(void* memory);
                MPI_Request mpi_request;
                void* memory;
            };
        public:
            enum direction { IN, OUT, LO };
            inline pipe(const packet_t& type, direction flow);
            inline ~pipe();
            inline size_t get_bound() const;
            inline request* create_request();
            inline request* attach_request(void* memory);
            inline void renew(request* r);
            inline void send(request* r);
            inline void recv(request* r);
            inline void spin();
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
        inline void ifetch(group* placement, size_t sid, size_t x, size_t y);

        inline group* world();
        inline void  init();
        inline void  finalize();
        inline std::pair<size_t*,size_t> id();
        inline packet_t& get_block_packet_type(size_t len);
        inline pipe* add_pipe(const packet_t& type, pipe::direction flow);
        inline pipe* get_pipe(const packet_t& type, pipe::direction flow);
        inline void  add_handler(const packet_t& type, void(*callback)(packet&));
        inline size_t get_volume() const;
        inline void  emit(packet* pack);
        inline void  spin();
        inline ~channel();

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
