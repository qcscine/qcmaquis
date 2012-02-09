#ifndef AMBIENT_GROUPS_MULTIRANK_H
#define AMBIENT_GROUPS_MULTIRANK_H
#include "ambient/channels/groups/group.h"
#include "ambient/utils/singleton.hpp"
#include <map>
#include <string>

#define UNDEFINED_RANK MPI_UNDEFINED

namespace ambient { namespace channels {

    class multirank : public singleton< multirank >
    {
    public:
        multirank();
    private: 
        std::map<std::string,int> map;
    public:
        int operator()(const group* grp) const;
        int operator()(const char* name = "ambient") const;
    };

} }

namespace ambient {
    extern channels::multirank& rank;
}
#endif
