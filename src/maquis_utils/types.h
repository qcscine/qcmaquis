
#ifndef MAQUIS_TYPES_TRAITS_H
#define MAQUIS_TYPES_TRAITS_H

namespace utils {
    template<class T> struct real_type { typedef T type; };
    template<class T> struct real_type<std::complex<T> > { typedef T type; };
}

#endif
