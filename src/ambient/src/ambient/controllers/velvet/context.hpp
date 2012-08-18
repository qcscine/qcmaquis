#include "ambient/channels/mpi/groups/group.h"
#include "ambient/utils/timings.hpp"

namespace ambient { namespace controllers { namespace velvet {


    inline context::context() : grp(NULL) { }

    inline context::~context(){
    }

    inline void context::set_group(group* grp){
        //this->functor->set_group(grp);
        this->grp = grp;
    }

} } }
