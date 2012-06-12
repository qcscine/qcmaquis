#ifndef AMBIENT_MODELS_VELVET_MEMSPEC
#define AMBIENT_MODELS_VELVET_MEMSPEC

namespace ambient { namespace models { namespace velvet {

    class memspec {
    public:
        inline memspec(size_t size);
        inline void* alloc() const;
        inline size_t get_bound() const;
        inline void* get_packet_t() const;
        size_t size;
    };

} } }

#endif
