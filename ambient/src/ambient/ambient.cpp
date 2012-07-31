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
// }}}
#include "ambient/utils/memory.h"
#include "ambient/channels/mpi/channel.h"
#include "ambient/models/velvet/model.h"
#include "ambient/controllers/velvet/controller.h"

namespace ambient {
    utils::bulk_memory& bulk_pool = utils::bulk_memory::instance();
    channels::mpi::multirank& rank = channels::mpi::multirank::instance();
    models::velvet::model& model = models::velvet::model::instance();
    channels::mpi::channel& channel = channels::mpi::channel::instance();
    controllers::velvet::controller& controller = controllers::velvet::controller::instance();
    controllers::velvet::context& ctxt = controllers::velvet::context::instance();
    io cout;
    io cerr;
}
