#ifndef AMBIENT_CONTROLLERS_VELVET_CONTEXT
#define AMBIENT_CONTROLLERS_VELVET_CONTEXT
#include "ambient/models/velvet/model.h"
#include "ambient/utils/singleton.hpp"

namespace ambient { namespace controllers { namespace velvet { 

    using ambient::channels::mpi::group;

    class context : public singleton< context > 
    {
    public:
        context();
       ~context();
    public:
// proxy functionality //
        context& operator()(int rank){ return *this;   } // proxy
        group* get_group()           { return grp;     }
        void set_group(group* grp);
    private:
        group* grp;
        cfunctor* functor;
// group class method duplicates
    public:
        int np,nq; //mask of the two cyclic distribution
        int  get_master()          { return ambient::rank.translate(grp->master, grp);   }
        bool involved()            { return ambient::rank.belongs(grp);                  }
        bool is_master()           { return ambient::rank.masters(grp);                  }
        int  get_rank()            { return grp->rank;                                   }
        int  get_size()            { return grp->size;                                   }
    };

} } }

namespace ambient {
    extern controllers::velvet::context& ctxt;
}

#endif
