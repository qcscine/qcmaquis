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
        inline void add_derivative(revision* r)        { r->set_generator(this); this->derivatives.push_back(r); }
        inline void add_dependency(revision* r)        { this->dependencies.push_back(r); }
        inline std::list<revision*>& get_derivatives() { return derivatives;  }
        inline std::list<revision*>& get_dependencies(){ return dependencies; }
        inline void set_group(group* g)                { grp = g; place();    }
        inline group* get_group()                      { return grp;          }
        void*  arguments[SFUNCTOR_ARITY];
        size_t revisions[SFUNCTOR_ARITY];
        std::list<revision*> derivatives;
        std::list<revision*> dependencies;
        group* grp;
    };

} } }

#endif
