#ifndef AMBIENT_CONTROLLERS_VELVET_CONTROLLER
#define AMBIENT_CONTROLLERS_VELVET_CONTROLLER

#include "ambient/controllers/velvet/cfunctor.h"
#include "ambient/controllers/velvet/iteratable.h"
#include "ambient/utils/collector.h"

#include <cilk/cilk.h>

namespace ambient { namespace controllers { namespace velvet {

    using ambient::models::velvet::revision;

    class controller : public singleton< controller >
    {
    public:
        controller();
       ~controller();
        void flush();
        void clear();
        void submit(cfunctor* f);
        void calloc(revision& r);
        void alloc (revision& r);
        void free  (revision& r);
        void sync  (revision& r);
        void sync  (revision& r, size_t target);
        template<typename T> void destroy(T* o);
    private:
        std::vector< cfunctor* > stack_m;
        std::vector< cfunctor* > stack_s;
        std::vector< cfunctor* >* chains;
        std::vector< cfunctor* >* mirror;
        collector garbage;
    };
    
} } }

namespace ambient {
    extern controllers::velvet::controller& controller;
}

#include "ambient/controllers/velvet/controller.hpp"
#include "ambient/utils/collector.hpp"
#include "ambient/controllers/velvet/iteratable.hpp"
#include "ambient/controllers/velvet/cfunctor.hpp"
#endif
