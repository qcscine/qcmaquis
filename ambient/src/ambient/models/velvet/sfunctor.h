#ifndef AMBIENT_MODELS_VELVET_SFUNCTOR
#define AMBIENT_MODELS_VELVET_SFUNCTOR
#include "ambient/channels/mpi/groups/group.h"

#define SFUNCTOR_ARITY 11

namespace ambient { namespace models { namespace velvet {
// spatial functor (modifier of objects)

    using ambient::channels::mpi::group;

    class sfunctor {
    public:
        virtual void deploy(size_t) = 0;
        void* arguments[SFUNCTOR_ARITY];
    };

} } }

#endif
