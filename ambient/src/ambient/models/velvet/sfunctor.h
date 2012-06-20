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
        virtual ~sfunctor(){
            this->grp->idle(); 
        } 
        //inline revision& get_vellum()       { return *vellum;   }
        //inline void set_vellum(revision& v) { vellum = &v;      }
        inline void set_group(group* g)     { grp = g; place(); }
        inline group* get_group()           { return grp;       }
        group* grp;
        //revision* vellum;
        void*  arguments[SFUNCTOR_ARITY];
        size_t revisions[SFUNCTOR_ARITY];
    };

} } }

#endif
