#ifndef AMBIENT_GROUPS_PACKET_MANAGER_H
#define AMBIENT_GROUPS_PACKET_MANAGER_H

#include "ambient/packets/types.h"
#include "ambient/packets/packet.h"
#include "ambient/auxiliary.h"
#include <list>

#define PULL_RESERVATION 20

using namespace ambient::packets; 

namespace ambient{ namespace groups{

    class group;
    class packet_manager
    {
    public: 
        enum locking_fsm { OPEN, LOOSE, CLOSURE, CLOSED };
        enum direction   { IN, OUT, LO };
    private: 
        packet_manager(packet_manager const&);             // copy constructor is private
        packet_manager& operator=(packet_manager const&);  // assignment operator is private
    public:
        packet_manager(groups::group* grp);
        class request
        {
        public:
            request(void* memory);
            MPI_Request mpi_request;
            void* memory;
            int fail_count;
        };
        class typed_q
        {
        public:
            typed_q(packet_manager* manager, const packet_t& type, direction flow, int reservation);
           ~typed_q();
            void spin();
            packet* get_target_packet();
            size_t get_active();
            request* attach_request(void* memory);
            request* get_request();
            void recv(request* r);
            void send(request* r);
            void delay(request* r);
        private:
            packet* target_packet;
            int reservation;
            std::list<request*> reqs;
        public:
            delegate packet_delivered;
            direction flow;
            packet_manager* manager;
            const packet_t& type;
        };
        typed_q* add_typed_q(const packet_t& type, direction flow, int reservation = 1);

        bool     subscribed(const packet_t& type);
        void     subscribe(const packet_t& type);
        void     add_handler(const packet_t& type, core::operation* callback);
        void     emit(packet* pack);
        typed_q* get_pipe(const packet_t& type, direction flow);

        void spin_loop();
        void spin();
        bool process_locking();
        groups::group* get_group();

        locking_fsm state;
        int closure_mutex;
        int approve_closure_mutex;
        std::list<typed_q*> qs;
    private:
        groups::group* grp;
        MPI_Comm* comm;
    };

} }
#endif
