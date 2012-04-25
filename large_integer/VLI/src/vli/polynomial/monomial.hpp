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

namespace vli {
    template <class Coeff>
    struct monomial {
        typedef unsigned int exponent_type;
        typedef Coeff coeff_type;
        typedef coeff_type value_type;

        /**
         * Constructor: Creates a monomial 1*J^j_exp*h^h_exp
         */
        explicit monomial(exponent_type j_exp = 0, exponent_type h_exp = 0)
        :j_exp_(j_exp), h_exp_(h_exp), coeff_(1) {
        }

        /*
        * These two functions are only called inside the vli_polynomial_gpu_function_hooks.hpp
        */
        value_type* p() {
            return coeff_.p();
        }
        
        value_type* p() const {
            return coeff_.p();
        }

        monomial operator - () const {
            monomial r(*this);
            r.coeff_ = -this->coeff_;
            return r;
        }
        monomial& operator *= (coeff_type const& c) {
            coeff_ *= c;
            return (*this);
        }
        
        monomial& operator *= (int c) {
            coeff_ *= c;
            return (*this);
        }
    
        void print(std::ostream& os) const {
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
    
    template <class Coeff>
    std::ostream& operator<<(std::ostream& os, monomial<Coeff> const& m){
        m.print(os);
        return os;
    }
    template <class Coeff, class T>
    monomial<Coeff> operator * (monomial<Coeff> m, T const& c) {
        m*=c;
        return m;
    }

    template <class Coeff, class T>
    monomial<Coeff> operator * (T const& c, monomial<Coeff> const& m) {
        return m*c;
    }

} //end namespace vli


#endif //VLI_MONOMIAL_HPP
