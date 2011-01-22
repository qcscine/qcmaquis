/*****************************************************************************
 *
 * ALPS DMFT Project - BLAS Compatibility headers
 *  Square Matrix Class
 *
 * Copyright (C) 2005 - 2010 by 
 *                              Emanuel Gull <gull@phys.columbia.edu>,
 *                              Brigitte Surer <surerb@phys.ethz.ch>
 *
 *
* This software is part of the ALPS Applications, published under the ALPS
* Application License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Application License along with
* the ALPS Applications; see the file LICENSE.txt. If not, the license is also
* available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/
 
#ifndef __ALPS_VECTOR_HPP__
#define __ALPS_VECTOR_HPP__

#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <alps/numeric/scalar_product.hpp>


namespace blas{
  template<typename T>
  class vector : public std::vector<T>
  {
    public:
      vector(std::size_t size=0, T const& initial_value = T())
      : std::vector<T>(size, initial_value)
      {
      }

      vector(vector const& v)
      : std::vector<T>(v)
      {
      }

      friend void swap(vector& x,vector& y)
      {
          std::swap(x, y);
      }
      
      inline T &operator()(const std::size_t i)
      {
          assert((i < this->size()));
          return this->operator[](i);
      }
      
      inline const T &operator()(std::size_t i) const 
      {
          assert((i < this->size()));
          return this->operator[](i);
      }
    
      vector& operator+=(const vector& rhs) 
      {
          assert(rhs.size() == this->size());
          plus_assign(this->begin(), this->end(), rhs.begin());
          return *this;
      }
      
      vector& operator-=(const vector& rhs) 
      {
          assert(rhs.size() == this->size());
          minus_assign(this->begin(), this->end(), rhs.begin());
          return *this;
      }
      
      vector& operator*=(T const& lambda) 
      {
          multiplies_assign(this->begin(), this->end(), lambda);
          return *this;
      }
  };  
    
    template<typename T>
    void insert(vector<T>& v, T value, std::size_t i)
    {
        assert((i <= v.size()));
        v.insert(v.begin()+i,value);
    }

    template <class InputIterator1, class InputIterator2>
    void plus_assign(InputIterator1 first1, InputIterator1 last1, InputIterator2 first2) 
    {
        std::transform(first1, last1, first2, first1, std::plus<typename std::iterator_traits<InputIterator2>::value_type >());
    }
   
    template<typename T>
    vector<T> operator+(vector<T> v1, const vector<T>& v2)
    {
        assert(v1.size() == v2.size());
        v1 += v2;
        return v1;
    }
    
    template <class InputIterator1, class InputIterator2>
    void minus_assign(InputIterator1 first1, InputIterator1 last1, InputIterator2 first2) 
    {
        std::transform(first1, last1, first2, first1, std::minus<typename std::iterator_traits<InputIterator2>::value_type >());
    }
    
    template<typename T>
    vector<T> operator-(vector<T> v1, const vector<T>& v2)  
    {
        assert(v1.size() == v2.size());
        v1 -= v2;
        return v1;
    }  

    template <class ForwardIterator, typename T>
    void multiplies_assign(ForwardIterator start1, ForwardIterator end1, T lambda) 
    {
        std::transform(start1, end1, start1, std::bind2nd(std::multiplies<T>(), lambda));
    }

    template <typename T>
    vector<T> operator * (T const& t, vector<T> v)
    {
        return v *= t;
    }
    
    template <typename T>
    vector<T> operator * (vector<T> v, T const& t)
    {
        return v *= t;
    }
    
    template <typename T>
    inline T scalar_product(const vector<T>& v1, const vector<T>& v2)
    {   
        return alps::numeric::scalar_product(v1,v2);
    }
    
    template <typename T>
    inline vector<T> exp(T c, vector<T> v)
    {
        using std::exp;
        vector<T> result(v.size());
        v*=c;
        std::transform(v.begin(), v.end(), result.begin(), static_cast<T(*)(T)> (&exp));
        return result;
    }
   
    inline vector<double> exp(double c, vector<double> v)
    {
        int s=v.size();
        vector<double> result(s);
        v*=c;
#ifdef VECLIB
        vecLib::vvexp(&result[0], &v[0], &s); 
#else
#ifdef ACML
        acml::vrda_exp(s, &v[0], &result[0]);
#else
#ifdef MKL
        mkl::vdExp(s,  &v[0], &result[0]);
#else
        using std::exp;
        std::transform(v.begin(), v.end(), result.begin(), static_cast<double(*)(double)> (&exp));
#endif
#endif
#endif  
        return result;
    }

  template <typename T>
  inline std::ostream &operator<<(std::ostream &os, const vector<T> &v)
  {
    os<<"[ ";
    for(unsigned int i=0;i<v.size()-1;++i){
      os<<v(i)<<", ";
    }
      os<< v(v.size()-1) << "]"<<std::endl;
    return os;
  }

    
} //namespace

#endif //__ALPS_VECTOR_HPP__
