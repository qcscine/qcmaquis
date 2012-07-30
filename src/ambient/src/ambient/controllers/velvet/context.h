#ifndef AMBIENT_CONTROLLERS_CONTEXT
#define AMBIENT_CONTROLLERS_CONTEXT
#include "ambient/models/velvet/model.h"
#include "ambient/utils/singleton.hpp"

namespace ambient { namespace controllers { namespace velvet { 

    using ambient::channels::mpi::group;
    using ambient::controllers::velvet::cfunctor;

    class context : public singleton< context > 
    {
    public:
        inline  context();
        inline ~context();
    public:
// proxy functionality //
        inline context& operator()(int rank){ return *this;   } // proxy
        inline group* get_group()           { return grp;     }
        inline void set_group(group* grp);
    private:
        group* grp;
        cfunctor* functor;
// group class method duplicates
    public:
        int np,nq; //mask of the two cyclic distribution
        inline int get_master_g()        { return grp->get_master_g();                         }
        inline int get_rank()            { return grp->get_rank();                             }
        inline int get_size()            { return grp->get_size();                             }
        inline const char* get_name()    { return grp->get_name();                             }
        inline bool involved()           { return grp->involved();                             }
        inline bool is_master()          { return grp->is_master();                            }
    };

} } }

namespace ambient {
    extern controllers::velvet::context& ctxt;
}

#endif
