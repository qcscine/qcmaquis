/*
 * Ambient Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *
 * Permission is hereby granted, free of charge, to any person or organization
 * obtaining a copy of the software and accompanying documentation covered by
 * this license (the "Software") to use, reproduce, display, distribute,
 * execute, and transmit the Software, and to prepare derivative works of the
 * Software, and to permit third-parties to whom the Software is furnished to
 * do so, all subject to the following:
 *
 * The copyright notices in the Software and this entire statement, including
 * the above license grant, this restriction and the following disclaimer,
 * must be included in all copies of the Software, in whole or in part, and
 * all derivative works of the Software, unless such copies or derivative
 * works are solely in the form of machine-executable object code generated by
 * a source language processor.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

namespace ambient { namespace models { namespace ssm {

    using ambient::memory::fixed;

    template<typename T> constexpr T op_sqrt (T a) { return std::sqrt(a); }
    template<typename T> constexpr T op_plus (T a, T b){ return a + b; }
    template<typename T> constexpr T op_minus(T a, T b){ return a - b; }
    template<typename T> constexpr T op_mul  (T a, T b){ return a * b; }
    template<typename T> constexpr T op_div  (T a, T b){ return a / b; }

    inline transformable::numeric_union::operator bool& (){ return b; }
    inline transformable::numeric_union::operator double& (){ return d; }
    inline transformable::numeric_union::operator std::complex<double>& (){ return c; }
    inline void transformable::numeric_union::operator = (bool value) { b = value; }
    inline void transformable::numeric_union::operator = (double value) { d = value; }
    inline void transformable::numeric_union::operator = (std::complex<double> value) { c = value; }

    inline void* transformable::operator new (size_t size, void* placement){ return placement; }
    inline void transformable::operator delete (void*, void*){ /* doesn't throw */ }

    template<typename T>
    inline transformable_value<T>::transformable_value(T value){
        this->v = value;
    }

    template<typename T>
    inline transformable::numeric_union transformable_value<T>::eval() const { 
        return this->v; 
    }

    template<typename T>
    inline transformable& transformable_value<T>::operator += (transformable& r){ 
        return *new(this)transformable_expr<T, decltype(&op_plus<T>), op_plus>(
                   &r, 
                   new (ambient::pool::calloc<fixed,sizeof_transformable()>()) transformable_value(*this)
               );
    }

    template <typename T, typename FP, FP OP>
    inline transformable::numeric_union transformable_expr<T,FP,OP>::eval() const {
        this->v = OP(this->l->eval(), this->r->eval());
        ambient::pool::free<fixed,sizeof_transformable()>((void*)this->l);
        ambient::pool::free<fixed,sizeof_transformable()>((void*)this->r);
        new((void*)this)transformable_value<T>(*(transformable_value<T>*)this);
        return this->v;
    }

    template <typename T, typename FP, FP OP>
    inline transformable& transformable_expr<T,FP,OP>::operator += (transformable& r){ 
        return *new(this)transformable_expr<T, decltype(&op_plus<T>), op_plus>(
                   &r, 
                   new (ambient::pool::calloc<fixed,sizeof_transformable()>()) transformable_expr<T,FP,OP>(*this)
               ); 
    }

    template <typename T, typename FP, FP OP>
    inline transformable_expr<T,FP,OP>::transformable_expr(const transformable* l){ 
        this->l = l; 
    }

    template <typename T, typename FP, FP OP>
    inline transformable_expr<T,FP,OP>::transformable_expr(const transformable* l, const transformable* r){ 
        this->l = l; 
        this->r = r; 
    }

} } }

