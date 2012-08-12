#ifndef AMBIENT_MODELS_VELVET_MODEL
#define AMBIENT_MODELS_VELVET_MODEL
#include "ambient/utils/singleton.hpp"
#include "ambient/utils/hashmap.hpp"
#include "ambient/utils/dim2.h"

#include "ambient/models/velvet/memspec.h"
#include "ambient/models/velvet/revision.h"
#include "ambient/models/velvet/history.h"
#include "ambient/models/velvet/sfunctor.h"

namespace ambient { namespace models { namespace velvet {

    class model : public singleton< model > {
    public:
        inline size_t time(const history* o);
        inline revision* add_revision(history* o, bool init = false);
        inline void insert(revision* r);
        inline revision* get_revision(size_t id) const;
    private:
        hashmap map;
    };

} } }

namespace ambient {
    extern models::velvet::model& model;
}

#include "ambient/models/velvet/model.hpp"
#include "ambient/models/velvet/history.hpp"
#include "ambient/models/velvet/revision.hpp"
#include "ambient/models/velvet/memspec.hpp"
#endif
