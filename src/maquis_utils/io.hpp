#ifndef MAQUIS_IO_HPP
#define MAQUIS_IO_HPP

namespace maquis {

#ifdef AMBIENT_IO
    extern ambient::io cout;
#else
    std::ostream& cout = std::cout;
#endif

}

#endif
