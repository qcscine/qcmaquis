#ifndef __MAQUIS_TYPES_TRAITS_HPP__
#define __MAQUIS_TYPES_TRAITS_HPP__

namespace maquis { namespace traits {

    template<class T> struct scalar_type { typedef typename T::value_type type; };
    template<class T> struct real_type { typedef typename real_type<typename T::value_type>::type type; };
    template<>        struct real_type<double> { typedef double type; };
    template<class T> struct real_type<std::complex<T> > { typedef T type; };
    template<class T> struct real_identity { static const T value; };
    template<class T> struct imag_identity { static const T value; };
    template<class T> struct real_identity<std::complex<T> > { static const std::complex<T> value; };
    template<class T> struct imag_identity<std::complex<T> > { static const std::complex<T> value; };
    template<class T> const T real_identity<T>::value = 1;
    template<class T> const T imag_identity<T>::value = 1;
    template<class T> const std::complex<T> real_identity<std::complex<T> >::value = std::complex<T>(1,0);
    template<class T> const std::complex<T> imag_identity<std::complex<T> >::value = std::complex<T>(0,1);

    template <class Matrix> struct transpose_view { typedef Matrix type; };
    template <class Matrix> struct transpose { typedef typename transpose_view<Matrix>::type type; };
    template <class Matrix> struct transpose< typename transpose_view<Matrix>::type > { typedef Matrix type; };

} }

#endif
