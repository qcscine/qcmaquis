#ifndef AMBIENT_CORE_WORKGROUP_CONTEXT_H
#define AMBIENT_CORE_WORKGROUP_CONTEXT_H
#include <utility>
#include <list>
#include <vector>

namespace ambient {

    class p_profile;

    class workgroup_context {
    public:
        workgroup_context();
        void bind(p_profile* pin, int i, int j = 0); // binds workgroup to the pin profile
        void discharge(core::operation* kernel);
    private:
        p_profile* pin;
        int i,j;
    };
}
#endif
