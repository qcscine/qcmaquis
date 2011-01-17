#include "ambient/groups/multirank.h"

namespace ambient{ namespace groups {

    multirank& multirank::instance(){
        static multirank* singleton = NULL;
        if(!singleton) singleton = new multirank();
        return *singleton;
    }
    multirank::multirank(){};

    void multirank::set(const group* grp, int rank){
        if(map.find(grp) == map.end()) map.insert(std::pair<const group*,int>(grp,rank));
        else map.find(grp)->second = rank;
    }

    int multirank::operator()(const group* grp){
        if(map.find(grp) == map.end()){ printf("Couldn't find requested group.\n"); return MPI_UNDEFINED; }
        else return map.find(grp)->second;
    }

    int multirank::operator()(const char* grp){
        if(map.find(group::group_map(grp)) == map.end()){ printf("Couldn't find requested group.\n"); return MPI_UNDEFINED; }
        else return map.find(group::group_map(grp))->second;
    }

} }
