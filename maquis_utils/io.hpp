#ifndef MAQUIS_IO_HPP
#define MAQUIS_IO_HPP

namespace maquis {

#ifdef AMBIENT_IO
    extern ambient::io cout;
    extern ambient::io cerr;
#else
    std::ostream& cout = std::cout;
    std::ostream& cerr = std::cerr;
#endif

}

#endif
