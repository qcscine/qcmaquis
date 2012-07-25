#ifndef AMBIENT_CONTROLLERS_CONTEXT
#define AMBIENT_CONTROLLERS_CONTEXT
#include "ambient/models/velvet/model.h"
#include "ambient/utils/singleton.hpp"

#define GET_TID ctxt.get_tid()
extern pthread_key_t pthread_tid;

namespace ambient { namespace controllers {     

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
        template<typename T>
        inline size_t get_revision_base(const T*);
        template<typename T>
        inline void set_revision_base(T* o, size_t base){
            o->set_thread_revision_base(base);
        }
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
        inline size_t get_tid()          { return *(size_t*)pthread_getspecific(pthread_tid);  }
        inline bool involved()           { return grp->involved();                             }
        inline bool is_master()          { return grp->is_master();                            }
        inline void set_tid(size_t);
    };

} }

namespace ambient {
    extern controllers::context& ctxt;
}

#endif
