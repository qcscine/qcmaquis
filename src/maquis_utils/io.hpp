#ifndef MAQUIS_IO_HPP
#define MAQUIS_IO_HPP

namespace maquis {

#ifdef AMBIENT_IO
    extern ambient::io cout;
    extern ambient::io cerr;
#else
    using std::cout;
    using std::cerr;
#endif

}

#endif
