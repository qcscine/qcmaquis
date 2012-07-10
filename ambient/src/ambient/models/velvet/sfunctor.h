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
        inline void set_derivative(revision* r) { derivative = r;     }
        inline revision* get_derivative()       { return derivative;  }
        inline void set_group(group* g)         { grp = g; place();   }
        inline group* get_group()               { return grp;         }
        void*  arguments[SFUNCTOR_ARITY];
        size_t revisions[SFUNCTOR_ARITY];
        revision* derivative;
        group* grp;
    };

} } }

#endif
