#ifndef MAQUIS_IO_HPP
#define MAQUIS_IO_HPP

#ifdef USE_AMBIENT
#include "ambient/ambient.hpp"
#endif

#include <iostream>
#include <string>
#include <iterator>

namespace maquis {
    
#ifdef USE_AMBIENT
    using ambient::cout;
    using ambient::cerr;
#else
    using std::cout;
    using std::cerr;
#endif
    
    template<class T, class Stream>
    class ostream_iterator_
    : public std::iterator<std::output_iterator_tag, void, void, void, void>
    {
    public:
        typedef std::string         char_type;
        typedef Stream              ostream_type;
        
    private:
        ostream_type*	    _M_stream;
        const char_type*	_M_string;
        
    public:
        ostream_iterator_(Stream& __s) : _M_stream(&__s), _M_string(0) {}
        
        ostream_iterator_(Stream& __s, const char_type& __c)
        : _M_stream(&__s), _M_string(&__c)  { }
        //
        // ostream_iterator(const ostream_iterator& __obj)
        // : _M_stream(__obj._M_stream), _M_string(__obj._M_string)  { }
        
        ostream_iterator_&
        operator=(const T& __value)
        {
            *_M_stream << __value;
            if (_M_string) *_M_stream << *_M_string;
            return *this;
        }
        
        ostream_iterator_&
        operator*()
        { return *this; }
        
        ostream_iterator_&
        operator++()
        { return *this; }
        
        ostream_iterator_&
        operator++(int)
        { return *this; }
    };
    
    template<class T, class Stream>
    ostream_iterator_<T, Stream> ostream_iterator(Stream& s) { return ostream_iterator_<T, Stream>(s); }
    template<class T, class Stream>
    ostream_iterator_<T, Stream> ostream_iterator(Stream& s, std::string const& c) { return ostream_iterator_<T, Stream>(s, c); }
    
}

#endif
