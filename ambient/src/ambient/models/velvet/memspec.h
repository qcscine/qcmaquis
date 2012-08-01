#ifndef AMBIENT_MODELS_VELVET_MEMSPEC
#define AMBIENT_MODELS_VELVET_MEMSPEC

namespace ambient { namespace models { namespace velvet {

    struct memspec {
        inline memspec(dim2 dim, size_t ts);
        inline void* alloc() const;
        inline void* calloc() const;
        inline size_t get_bound() const;
        inline void* get_packet_t() const;
        dim2 dim;
        size_t size;
    };

} } }

#endif
