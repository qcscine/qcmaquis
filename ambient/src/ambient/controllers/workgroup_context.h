#ifndef AMBIENT_CORE_WORKGROUP_CONTEXT_H
#define AMBIENT_CORE_WORKGROUP_CONTEXT_H
#include <utility>
#include <list>
#include <vector>

namespace ambient { namespace controllers {

    class workgroup_context {
    public:
        workgroup_context();
        void bind(models::operation* op); // binds workgroup to the pin object
        void discharge(models::operation* kernel);
    private:
        models::imodel::object* pin;
        models::operation* bind_op;
    };

} }
#endif
