#ifndef AMBIENT_MODELS_VELVET_SFUNCTOR
#define AMBIENT_MODELS_VELVET_SFUNCTOR
#include "ambient/channels/mpi/groups/group.h"

#define SFUNCTOR_ARITY 10

namespace ambient { namespace models { namespace velvet {

    using ambient::channels::mpi::group;

    // spatial functor (modifier of objects)
    class sfunctor {
    public:
        virtual void place() = 0;
        void set_group(group* g)                  { grp = g; place();    }
        group* get_group()                        { return grp;          }
        void*  arguments[SFUNCTOR_ARITY];
        group* grp;
    };

} } }

#endif
