/*
 *  monome.h
 *
 *  Created by Tim Ewart (timothee.ewart@unige.ch) and  Andreas Hehn (hehn@phys.ethz.ch) on 18.03.11.
 *  Copyright 2011 University of Geneva and Eidgenössische Technische Hochschule Züric. All rights reserved.
 *
 */

#ifndef VLI_MONOMIAL_HPP
#define VLI_MONOMIAL_HPP

#include <ostream>
#include <cmath>

namespace vli
{
    template <class Vli>
    struct monomial
    {
        typedef unsigned int exponent_type;
        typedef Vli value_type;

        /**
         * Constructor: Creates a monomial 1*J^j_exp*h^h_exp
         */
        explicit monomial(exponent_type j_exp = 0, exponent_type h_exp = 0)
        :j_exp_(j_exp), h_exp_(h_exp), coeff_(1){
        }

        explicit monomial(Vli a, exponent_type j_exp = 0, exponent_type h_exp = 0)
        :j_exp_(j_exp), h_exp_(h_exp), coeff_(a){
        }
        
        /*
        * These two functions are only called inside the vli_polynomial_gpu_function_hooks.hpp
        */
//        value_type* p(){
//            return coeff_.p();
//        }
//        
//        value_type* const p() const{
//            return coeff_.p();
//        }

        monomial& operator *= (Vli const& c){
            coeff_ *= c;
            return (*this);
        }
        
        monomial& operator *= (int c){
            coeff_ *= c;
            return (*this);
        }

        /**
        * For the GPU proxy 
        * MEGA DANGEROUS
        * Andreas one free drink, if you get a safer solution solution
        * I tried to specify but my compilo is blind
        */        
//      template<class Order>
//      monomial& operator *= (typename vector_polynomial_gpu<polynomial_gpu<Vli,Order> >::proxy const& p){
 /*
        template <class T>                                    
        monomial& operator *= (T const& p){
            Vli Tmp = p.BuildProxyToVli();
            coeff_ *= Tmp;
            return (*this);
        }
       */     
        void print(std::ostream& os) const
        {
            if(coeff_ > 0)
                os<<"+";
            os<<coeff_<<"*J^"<<j_exp_<<"*h^"<<h_exp_;
        }

        bool operator == (monomial const& m) const{
            return (j_exp_ == m.j_exp_) && (h_exp_ == m.h_exp_) && (coeff_ == m.coeff_);
        }
    
        exponent_type j_exp_;
        exponent_type h_exp_;
        value_type coeff_;
    };
    
    template<class Vli>
    std::ostream& operator<<(std::ostream& os, monomial<Vli> const& m){
        m.print(os);
        return os;
    }
    template <typename Vli, typename T>
    monomial<Vli> operator * (monomial<Vli> m, T const& c)
    {
        m*=c;
        return m;
    }

    template <typename Vli, typename T>
    monomial<Vli> operator * (T const& c, monomial<Vli> const& m)
    {
        return m*c;
    }

}


#endif //VLI_MONOMIAL_HPP
