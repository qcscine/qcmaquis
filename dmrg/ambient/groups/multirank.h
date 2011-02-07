#ifndef AMBIENT_GROUPS_MULTIRANK_H
#define AMBIENT_GROUPS_MULTIRANK_H
#include <mpi.h>
#include "ambient/groups/group.h"

namespace ambient{ namespace groups{

    class multirank
    {
    private: 
        multirank();                               // constructor is private
        multirank(multirank const&){};             // copy constructor is private
        multirank& operator=(multirank const&){};  // assignment operator is private
        std::map<std::string,int> map;
    public:
        static multirank& instance();
    public:
        void set(const group* grp, int rank);
        int operator()(const group* grp);
        int operator()(const char* grp = "ambient");
    };

} }
#endif
