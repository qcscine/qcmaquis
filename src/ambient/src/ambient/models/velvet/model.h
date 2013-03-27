#ifndef AMBIENT_MODELS_VELVET_MODEL
#define AMBIENT_MODELS_VELVET_MODEL
#include "ambient/utils/singleton.hpp"
#include "ambient/utils/dim2.h"

#include "ambient/models/velvet/memspec.h"
#include "ambient/models/velvet/revision.h"
#include "ambient/models/velvet/history.h"
#include "ambient/models/velvet/transformable.h"

#define MAX_SID 2147483647

namespace ambient { namespace models { namespace velvet {

    class model : public singleton< model > {
    public:
        model() : clock(0), sid(0) {}
        template<ambient::rstate S> 
        void add_revision(history* o, void* g = NULL);
        void use_revision(history* o);
        bool feeds(const revision* r);
        bool common(const revision* r);
        size_t time(const history* o);
        void touch(const history* o);
        void index(revision* r);
        void index(transformable* v);
        size_t clock;
    private:
        int sid;
    };

} } }

namespace ambient {
    extern models::velvet::model& model;
}

#include "ambient/models/velvet/model.hpp"
#include "ambient/models/velvet/transformable.hpp"
#include "ambient/models/velvet/history.hpp"
#include "ambient/models/velvet/revision.hpp"
#include "ambient/models/velvet/memspec.hpp"
#endif
