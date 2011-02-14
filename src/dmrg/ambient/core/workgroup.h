#ifndef AMBIENT_CORE_WORKGROUP_H
#define AMBIENT_CORE_WORKGROUP_H
#include <utility>
#include <list>
#include <vector>

namespace ambient {

    class p_profile;

    class workgroup {
    public:
        workgroup(p_profile** p, int i, int j = 0, int k = 0);
        void* item(int i, int j = 0, int k = 0);
        p_profile** profile;
        void* data;
        bool initialized;
        int owner;
        int i, j, k;
    };
}
#endif
