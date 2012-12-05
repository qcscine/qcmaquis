#ifndef AMBIENT_MODELS_VELVET_MEMSPEC
#define AMBIENT_MODELS_VELVET_MEMSPEC

namespace ambient { namespace models { namespace velvet {

    struct memspec {
        memspec(dim2 dim, size_t ts);
        size_t size;
        dim2 dim;
    };

} } }

#endif
