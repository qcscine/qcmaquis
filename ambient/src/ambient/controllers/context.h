#ifndef AMBIENT_CONTROLLERS_CONTEXT
#define AMBIENT_CONTROLLERS_CONTEXT
#include "ambient/models/velvet/model.h"
#include "ambient/utils/singleton.hpp"

extern pthread_key_t pthread_tid;

#define GET_TID 0 // this->get_tid()

namespace ambient { namespace controllers {     

    using ambient::channels::mpi::group;
    using ambient::controllers::velvet::cfunctor;
    using ambient::controllers::velvet::iteratable;

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
        inline size_t get_revision_base(const iteratable<T>*);
        template<typename T>
        inline void set_revision_base(iteratable<T>* o, size_t base){
            o->set_thread_revision_base(base);
        }
    private:
        group* grp;
        cfunctor* functor;
        dim2 thread_block_id[AMBIENT_THREADS_LIMIT];
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
        inline void set_block_id(dim2 k) { thread_block_id[GET_TID] = k;                       }
        inline dim2 get_block_id()       { return thread_block_id[GET_TID];                    }
        inline void set_tid(size_t);
    };

} }

namespace ambient {
    extern controllers::context& ctxt;
}

#endif
