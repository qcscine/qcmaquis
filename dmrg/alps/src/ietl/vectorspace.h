/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2002 by Prakash Dayal <prakash@comp-phys.org>,
*                            Matthias Troyer <troyer@comp-phys.org>
*
* This software is part of the ALPS libraries, published under the ALPS
* Library License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Library License along with
* the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
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

/* $Id: vectorspace.h,v 1.13 2004/06/29 08:31:02 troyer Exp $ */

#ifndef IETL_VECTORSPACE__H
#define IETL_VECTORSPACE__H
#include <ietl/traits.h>
#include <boost/smart_ptr.hpp>

namespace ietl {
  template<class V>
    class vectorspace {
    public:
    typedef V vector_type;
    typedef typename V::value_type scalar_type;
    typedef typename V::size_type size_type;
    
    vectorspace(size_type n):n_(n){}
    
    inline size_type vec_dimension() const {
      return n_;
    }
    vector_type new_vector() const {
      vector_type v(n_);
      clear(v);
      return v;
    }
    
    void project(vector_type&) const {
    }
    
    private:
    size_type n_;
  };
  
  template <class V, class S> class scaled_vector_wrapper;
  
  template <class V>
    class vector_wrapper : public boost::shared_ptr<V> {
    typedef boost::shared_ptr<V> super_type;
    public:
    vector_wrapper(V* p) : boost::shared_ptr<V>(p) {}
    operator V& () { return *super_type::get();}
    operator const V& () const { return *super_type::get();}
    const vector_wrapper operator += (const vector_wrapper& x) { *super_type::get() += *x.get(); return *this;}
    const vector_wrapper operator -= (const vector_wrapper& x) { *super_type::get() -= *x.get(); return *this;}
    template <class T> const vector_wrapper& operator *= (T x) { *super_type::get() *= x; return *this;}
    template <class T> const vector_wrapper& operator /= (T x) { *super_type::get() /= x; return *this;}
    template <class S>
    const vector_wrapper& operator += (const scaled_vector_wrapper<V,S>& x) 
    { *super_type::get() += x.scalar()*x.vector(); return *this;}
    template <class S>
    const vector_wrapper& operator -= (const scaled_vector_wrapper<V,S>& x) 
    { *super_type::get() -= x.scalar()*x.vector(); return *this;}
    template <class S>
    const vector_wrapper& operator = (const scaled_vector_wrapper<V,S>& x) 
    { *super_type::get() = x.scalar()*x.vector(); return *this;}
      };
  
  template<class VS>
    void project(typename ietl::vectorspace_traits<VS>::vector_type& v, const VS& vs) {
    vs.project(v);
  }
  
  template<class V>
    class wrapper_vectorspace {
    public:
    typedef vector_wrapper<V> vector_type;
    typedef typename V::value_type scalar_type;
    typedef typename V::size_type size_type;
    
    wrapper_vectorspace(size_type n):n_(n){}
    
    inline size_type vec_dimension() const{
      return n_;
    }
    vector_type new_vector() const {
      vector_wrapper<V> v(new V(n_));
      clear(*v.get());
      return v;
    }
    
    void project(vector_type& src) const {
    }  
    private:
    size_type n_;
  };
  
template <class VS>
  typename ietl::vectorspace_traits<VS>::vector_type new_vector(const VS& vs) {
  return vs.new_vector();
}

 template <class VS>
   typename ietl::vectorspace_traits<VS>::size_type vec_dimension(const VS& vs) {
   return vs.vec_dimension();
 } 

   /*
 template <class V, class VS>
   void project(ietl::vector_wrapper<V>& v, const VS& vs) {
   vs.project(v);
 }
 */

 template <class V, class S>
   class scaled_vector_wrapper {
   public:
   scaled_vector_wrapper(const ietl::vector_wrapper<V>& v, S s)
     : v_(v), s_(s)
     {}
   
   const V& vector() const { return *v_.get();}
   S scalar() const { return s_;}
   private:
   const ietl::vector_wrapper<V>& v_;
   S s_;
 }; 
} // matches namespace ietl

// wrapper forwarders

namespace ietl {
  template <class V>
    void copy(const ietl::vector_wrapper<V>& src, ietl::vector_wrapper<V>& dst) {
    ietl::copy(*src.get(),*dst.get());
  } 
  

  template <class V>
  typename number_traits<typename V::value_type>::magnitude_type two_norm(const ietl::vector_wrapper<V>& src) {
    return ietl::two_norm(*src.get());
  }
  
 template <class V>
   typename V::value_type dot(const ietl::vector_wrapper<V>& src1, const ietl::vector_wrapper<V>& src2) {
   return ietl::dot(*src1.get(),*src2.get());
 }
 
 template <class A, class V>
   void mult(A a, const ietl::vector_wrapper<V>& src, ietl::vector_wrapper<V>& dst) {
   ietl::mult(a,*src.get(),*dst.get());
 } 

  template <class V, class GEN>
   void generate(ietl::vector_wrapper<V>& src, const GEN& gen) {
   ietl::generate(*src.get(),gen);
 }
 
 template <class V, class S>
   ietl::scaled_vector_wrapper<V,S> operator*(const ietl::vector_wrapper<V>& v, S s) {
   return ietl::scaled_vector_wrapper<V,S>(v,s);
 }

 template <class V, class S>
   ietl::scaled_vector_wrapper<V,S> operator*(S s, const ietl::vector_wrapper<V>& v) {
   return ietl::scaled_vector_wrapper<V,S>(v,s);
 }

 template <class V, class S>
   ietl::scaled_vector_wrapper<V,S> operator/(const ietl::vector_wrapper<V>& v, S s) {
   return ietl::scaled_vector_wrapper<V,S>(v,1./s);
 }
 
} // end of namespace ietl.
#endif
