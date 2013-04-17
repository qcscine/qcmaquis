#ifndef AMBIENT_NUMERIC_TRAITS
#define AMBIENT_NUMERIC_TRAITS

namespace ambient { namespace numeric { namespace traits {

    template<class T> struct real_type { typedef typename real_type<typename T::value_type>::type type; };
    template<>        struct real_type<double> { typedef double type; };
    template<class T> struct real_type<std::complex<T> > { typedef T type; };

} } }

#endif
