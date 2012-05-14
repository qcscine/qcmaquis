#ifndef AMBIENT_INTERFACE_CHANNEL
#define AMBIENT_INTERFACE_CHANNEL

namespace ambient { namespace channels {

    class group;
    class packet_t;
    class ichannel {
    public:
        class pipe {
            enum direction { IN, OUT, LO };
        };
        class packet {
        public:
            virtual const packet_t& get_t() = 0;
            virtual size_t get_bound(size_t field) = 0;
            virtual void* get_memory() = 0;
        };
        virtual size_t get_volume() const = 0;
        virtual void init() = 0;
        virtual void finalize() = 0;
        virtual void emit(packet* p) = 0;
        virtual void spin() = 0;
        virtual void ifetch(channels::group* placement, size_t gid, size_t sid, size_t i, size_t j) = 0;
        virtual std::pair<size_t*,size_t> id() = 0;
        virtual channels::packet_t& get_block_packet_type(size_t) = 0;
        virtual group* world() = 0;
    };
    
} }

namespace ambient {
    extern channels::ichannel& channel;
}

#endif
