#ifndef AMBIENT_INTERFACE
#define AMBIENT_INTERFACE
#ifndef AMBIENT
#define AMBIENT
#endif
#define AMBIENT_THREADS 6
#define AMBIENT_THREADS_LIMIT 12
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
#include <algorithm>
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
#include "ambient/interface/future.hpp"

namespace ambient{

    using models::velvet::history;
    using models::velvet::revision;
    using controllers::velvet::iteratable;
    using controllers::velvet::c_revision;
    using controllers::velvet::w_revision;
    using controllers::velvet::p_revision;
    using controllers::velvet::r_revision;
    using controllers::velvet::cfunctor;

    template <typename T> inline revision&   ui_l_current(T& obj){ return obj.ui_l_revision_0(); } // (logistics)
    template <typename T> inline c_revision& ui_c_current(T& obj){ return obj.ui_c_revision_0(); } // (computation)
    template <typename T> inline w_revision& ui_w_updated(T& obj){ return obj.ui_w_revision_1(); } // (reusable weak)
    template <typename T> inline p_revision& ui_p_updated(T& obj){ return obj.ui_p_revision_1(); } // (purged malloc)
    template <typename T> inline r_revision& ui_r_updated(T& obj){ return obj.ui_r_revision_1(); } // (reusable)

    template<typename T> inline dim2   ui_c_get_dim        (T& ref){ return ref.spec.dim;   }
    template<typename T> inline size_t ui_c_get_mem_size   (T& ref){ return ref.spec.size;  }
}

#endif
