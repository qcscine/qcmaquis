#ifndef REDI_INDEX_TUPLE_H
#define REDI_INDEX_TUPLE_H

// Copyright Jonathan Wakely 2012
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

namespace redi
{
  /// A type that represents a parameter pack of zero or more integers.
  template<unsigned... Indices>
    struct index_tuple
    {
      /// Generate an index_tuple with an additional element.
      template<unsigned N>
        using append = index_tuple<Indices..., N>;
    };

  /// Unary metafunction that generates an index_tuple containing [0, Size)
  template<unsigned Size>
    struct make_index_tuple
    {
      typedef typename make_index_tuple<Size-1>::type::template append<Size-1>
        type;
    };

  // Terminal case of the recursive metafunction.
  template<>
    struct make_index_tuple<0u>
    {
      typedef index_tuple<> type;
    };

  template<typename... Types>
    using to_index_tuple = typename make_index_tuple<sizeof...(Types)>::type;

}  // namespace redi

#endif  // REDI_INDEX_TUPLE_H

// vi: set ft=cpp sw=2:
