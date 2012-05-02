#ifndef MAQUIS_TYPES_TRAITS_H
#define MAQUIS_TYPES_TRAITS_H

namespace utils {

    template<class T> struct real_type { typedef typename real_type<typename T::value_type>::type type; };
    template<>        struct real_type<double> { typedef double type; };
    template<class T> struct real_type<std::complex<T> > { typedef T type; };
    
    template<class T> struct real_identity { static const T value; };
    template<class T> struct imag_identity { static const T value; };
    
    template<class T> struct real_identity<std::complex<T> > { static const std::complex<T> value; };
    template<class T> struct imag_identity<std::complex<T> > { static const std::complex<T> value; };
    
    // complex identities
    
    template<class T>
    const T real_identity<T>::value = 1;
    template<class T>
    const T imag_identity<T>::value = 1;
    
    template<class T>
    const std::complex<T> real_identity<std::complex<T> >::value = std::complex<T>(1,0);
    template<class T>
    const std::complex<T> imag_identity<std::complex<T> >::value = std::complex<T>(0,1);

}

#endif

