#ifndef AMBIENT_CORE_OPERATION_H
#define AMBIENT_CORE_OPERATION_H
#include <stdlib.h>
#define pinned ambient::core::ambient_pin* , 
#define marked NULL,

namespace ambient{
    class p_profile; 
    namespace groups { class group; }
}
namespace ambient{ namespace core{

    class ambient_pin{};

    class operation{
    public:
        #include "ambient/core/operation/operation.pp.h"
        void init();     // some init that doesn't depend upon arguments
        void perform();  // executes operation
        void invoke();   // executes operation (clean way)
        void set_scope(groups::group* scope);
        groups::group* get_scope();
        void preprocess();
        void finalize();
        void(operation::*prototype)();
        void(*operation_ptr)();
        groups::group* scope;
        void** arguments;
        p_profile** profiles;
        p_profile* pin;
        size_t count;
        bool is_extracted;
    };

    class out_of_scope_e{};

} }

#endif
