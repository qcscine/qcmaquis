#ifndef AMBIENT_MODELS_VELVET_MEMSPEC
#define AMBIENT_MODELS_VELVET_MEMSPEC

namespace ambient { namespace models { namespace velvet {

    struct memspec {
        memspec(dim2 dim, size_t ts);
        void* alloc() const;
        void* calloc() const;
        void free(void*) const;
        size_t get_bound() const;
        void* get_packet_t() const;
        size_t size;
        dim2 dim;
    };

} } }

#endif
