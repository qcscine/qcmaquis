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
        template<typename T> inline size_t time(const T* o);
        template<typename T> inline revision& add_revision(T* o);
        inline bool is_atomic(const history* o);
        inline size_t get_block_lda(history* o);
        inline dim2 get_current_dim(const history* o);
        inline void set_current_dim(history* o, dim2);
        inline void insert(layout* l);
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
