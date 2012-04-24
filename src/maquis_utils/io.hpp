#ifndef MAQUIS_IO_HPP
#define MAQUIS_IO_HPP

namespace maquis {

#ifdef AMBIENT
    extern ambient::cout cout;
#else
    extern std::cout cout;
#endif

}

#endif
