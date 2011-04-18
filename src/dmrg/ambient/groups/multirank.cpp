#include "ambient/mpi.h"
#include "ambient/ambient.h"
#include "ambient/groups/multirank.h"

namespace ambient{ namespace groups {

    multirank& multirank::instance(){
        static multirank* singleton = NULL;
        if(!singleton) singleton = new multirank();
        return *singleton;
    }
    multirank::multirank(){};

    int multirank::operator()(const group* grp) const {
        return grp->rank;
    }

    int multirank::operator()(const char* name) const {
        group* grp = groups::group_map(name);
        if(grp != NULL) return grp->rank;
        return MPI_UNDEFINED;
    }

    bool multirank::is_master(const char* name) const
    {
          group* grp = groups::group_map(name);
          return is_master(grp);
    }

    bool multirank::is_master(const group* grp) const
    {
          if((*this)(grp) == UNDEFINED_RANK) return false;
          return grp->translate_up_rank((*this)(grp)) == grp->translate_up_rank(grp->master);
    }

} }
