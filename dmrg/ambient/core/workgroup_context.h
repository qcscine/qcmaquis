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
        void bind(core::operation* op); // binds workgroup to the pin profile
        void discharge(core::operation* kernel);
        void finalize();
    private:
        p_profile* pin;
        core::operation* bind_op;
        int i,j;
    };
}
#endif
