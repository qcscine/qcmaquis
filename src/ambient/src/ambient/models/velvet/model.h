#ifndef AMBIENT_MODELS_VELVET_MODEL
#define AMBIENT_MODELS_VELVET_MODEL
#include "ambient/utils/singleton.hpp"
#include "ambient/utils/hashmap.hpp"
#include "ambient/utils/dim2.h"

#include "ambient/models/velvet/layout.h"
#include "ambient/models/velvet/sfunctor.h"
#include "ambient/models/velvet/revision.h"
#include "ambient/models/velvet/reduction.h"
#include "ambient/models/velvet/history.h"

namespace ambient { namespace models { namespace velvet {

    class model : public singleton< model > {
    public:
        inline model();
        inline void insert(layout* l);
        inline void add_revision(history* obj);
        inline void update_revision(revision* r, group* placement);
        inline layout* get_layout(size_t id) const;
        inline model& operator>>(dim2);
        inline model& operator, (dim2);
        inline ~model();
        dim2 mem_dim;
    private:
        hashmap map;
    };

} } }

namespace ambient {
    extern models::velvet::model& model;
}

#include "ambient/models/velvet/model.hpp"
#include "ambient/models/velvet/history.hpp"
#include "ambient/models/velvet/reduction.hpp"
#include "ambient/models/velvet/revision.hpp"
#include "ambient/models/velvet/layout.hpp"
#endif
