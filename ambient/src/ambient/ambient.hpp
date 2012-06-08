#ifndef AMBIENT_INTERFACE
#define AMBIENT_INTERFACE
#ifndef AMBIENT
#define AMBIENT
#endif
// {{{ system includes
#include <mpi.h>
#include <complex>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string>
#include <limits>
#include <vector>
#include <set>
#include <map>
#include <list>
#include <memory.h>
#include <stdarg.h>
#include <ctype.h>
#include <iostream>
#include <fstream>
#include <sys/time.h>
#include <pthread.h>
#define BOOST_SP_NO_SP_CONVERTIBLE
#include <boost/intrusive_ptr.hpp>
#include <boost/shared_ptr.hpp>
// }}}
#include "ambient/channels/mpi/channel.h"
#include "ambient/models/velvet/model.h"
#include "ambient/controllers/velvet/controller.h"
#include "ambient/utils/auxiliary.hpp"
#include "ambient/utils/memory.hpp"
#include "ambient/utils/io.hpp"
#include "ambient/interface/typed.hpp"
#include "ambient/interface/kernel.hpp"
#include "ambient/interface/parallel.hpp"
#include "ambient/interface/future.hpp"

namespace ambient{

    using models::velvet::history;
    using models::velvet::revision;
    using controllers::velvet::iteratable;
    using controllers::velvet::fast_revision;
    using controllers::velvet::slow_revision;
    using controllers::velvet::cfunctor;

    template <typename T> inline revision&      ui_l_current(T& obj){ return obj.ui_l_revision_0();     } // (logistics)
    template <typename T> inline revision&      ui_l_updated(T& obj){ return obj.ui_l_revision_1();     } // (logistics, debug)
    template <typename T> inline fast_revision& ui_c_current(T& obj){ return obj.ui_c_revision_0();     } // (computation)
    template <typename T> inline slow_revision& ui_c_updated(T& obj){ return obj.ui_c_revision_1();     } // (memcheck)

    inline void pin(cfunctor* o, revision& r, int x, int y){ r.content->entries[x][y]->assignments.push_back(o); }
    inline void assign(revision& r, int x, int y){ ambient::controller.ifetch_block(r, x, y); }

    template<typename T> inline dim2 ui_c_get_dim        (T& ref){ return ui_c_current(ref).get_layout().dim;            }
    template<typename T> inline dim2 ui_c_get_mem_dim    (T& ref){ return ui_c_current(ref).get_layout().mem_dim;        }
    template<typename T> inline dim2 ui_c_get_grid_dim   (T& ref){ return ui_c_current(ref).get_layout().grid_dim;       }
    template<typename T> inline dim2 ui_c_get_dim_u      (T& ref){ return ui_c_updated(ref).get_layout().dim;            }
    template<typename T> inline dim2 ui_c_get_mem_dim_u  (T& ref){ return ui_c_updated(ref).get_layout().mem_dim;        }
    template<typename T> inline dim2 ui_c_get_grid_dim_u (T& ref){ return ui_c_updated(ref).get_layout().grid_dim;       }

    inline dim2 ui_l_get_grid_dim(revision& r){ return r.content->grid_dim; }

}

#endif
