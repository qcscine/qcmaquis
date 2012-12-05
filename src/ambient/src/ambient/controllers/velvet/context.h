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
        int get_master_g()        { return grp->get_master_g();                         }
        int get_rank()            { return grp->get_rank();                             }
        int get_size()            { return grp->get_size();                             }
        const char* get_name()    { return grp->get_name();                             }
        bool involved()           { return grp->involved();                             }
        bool is_master()          { return grp->is_master();                            }
    };

} } }

namespace ambient {
    extern controllers::velvet::context& ctxt;
}

#endif
