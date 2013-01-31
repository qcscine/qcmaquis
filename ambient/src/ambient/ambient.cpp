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
#include "ambient/utils/memory.hpp"
#include "ambient/channels/mpi/channel.h"
#include "ambient/models/velvet/model.h"
#include "ambient/controllers/velvet/controller.h"

namespace ambient {
    memory::bulk& bulk = memory::bulk::instance();
    memory::pool& pool = memory::pool::instance();
    models::velvet::model& model = models::velvet::model::instance();
    channels::mpi::channel& channel = channels::mpi::channel::instance();
    channels::mpi::multirank& rank = channels::mpi::multirank::instance();
    controllers::velvet::controller& controller = controllers::velvet::controller::instance();
    utils::fstream fout;
    utils::mpostream cout;
    utils::mpostream cerr;
}
