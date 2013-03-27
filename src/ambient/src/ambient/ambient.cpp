#include "ambient/ambient.hpp"

namespace ambient {
    int scope<single>::factor = 1;
    int scope<single>::effect = 0;
    int scope<single>::iterator = 0;
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
