#include "ambient/ambient.h"
#include "ambient/groups/multirank.h"

namespace ambient{ namespace groups {

    multirank& multirank::instance(){
        static multirank* singleton = NULL;
        if(!singleton) singleton = new multirank();
        return *singleton;
    }
    multirank::multirank(){};

    void multirank::set(const group* grp, int rank){
        if(map.find(grp->name) == map.end()) map.insert(std::pair<std::string,int>(grp->name,rank));
        else map.find(grp->name)->second = rank;
    }

    int multirank::operator()(const group* grp) const {
        if(map.find(grp->name) == map.end()){ printf("Couldn't find requested group.\n"); return MPI_UNDEFINED; }
        else return map.find(grp->name)->second;
    }

    int multirank::operator()(const char* grp) const {
        if(map.find(grp) == map.end()){ printf("Couldn't find requested group.\n"); return MPI_UNDEFINED; }
        else return map.find(grp)->second;
    }

    bool multirank::is_master(const char* grp_name) const
    {
          groups::group* grp = groups::group_map(grp_name);
          return is_master(grp);
    }

    bool multirank::is_master(const group* grp) const
    {
          if((*this)(grp) == UNDEFINED_RANK) return false;
          return grp->translate_rank((*this)(grp)) == grp->translate_rank(grp->master);
    }

} }
