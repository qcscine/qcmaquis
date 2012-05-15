#ifndef AMBIENT_CONTROLLERS_CONTEXT
#define AMBIENT_CONTROLLERS_CONTEXT
#include "ambient/models/v_model.h"
#include "ambient/utils/singleton.hpp"

extern pthread_key_t pthread_env;
extern pthread_key_t pthread_tid;

namespace ambient { namespace controllers {

    class context : public singleton< context > 
    { // scalable multiprocessor
    public:
        context();
    public:
// proxy functionality //
        inline context& operator()(const int rank){ return *this; } // ??
        void set_group(channels::group* grp);
        inline channels::group* get_group(){ return this->grp; }
        inline void set_op(models::imodel::modifier* op){ this->op = op; }
        inline models::imodel::modifier* get_op(){ return this->op; }
        size_t get_revision_base(const models::imodel::object*);
        void set_revision_base(models::imodel::object*, size_t);
    private:
        channels::group* grp;
        models::imodel::modifier* op;
// proxy functionality //
// group class method duplicates
    public:
        enum { MARKUP, EXECUTE } state;
        int np,nq; //mask of the two cyclic distribution
        inline int get_master_g()    { return this->grp->get_master_g(); }
        inline int get_rank()        { return this->grp->get_rank();     }
        inline int get_size()        { return this->grp->get_size();     }
        inline const char* get_name(){ return this->grp->get_name();     }
        inline dim2 get_block_id()   { return *(dim2*)pthread_getspecific(pthread_env); }
        inline size_t get_tid()      { return *(size_t*)pthread_getspecific(pthread_tid); }
        inline bool involved()       { return (this->state == MARKUP ? false : this->grp->involved()); }
        inline bool is_master()      { return this->grp->is_master();    }
        void set_block_id(dim2);
        void set_tid(size_t);
    };

} }

namespace ambient {
    extern controllers::context& ctxt;
}

#endif
