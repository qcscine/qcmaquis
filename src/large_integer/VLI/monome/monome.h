/*
 *  monome.h
 *
 *  Created by Tim Ewart (timothee.ewart@unige.ch) and  Andreas Hehn (hehn@phys.ethz.ch) on 18.03.11.
 *  Copyright 2011 University of Geneva and Eidgenössische Technische Hochschule Züric. All rights reserved.
 *
 */

#ifndef VLI_MONOME_H
#define VLI_MONOME_H

#include <ostream>
#include <cmath>

namespace vli
{	
	template <class Vli>
    struct monomial
    {
        typedef typename Vli::size_type size_type;      // Type of the exponents (has 
         
        /**
         * Constructor: Creates a monomial 1*J^j_exp*h^h_exp
         */
        explicit monomial(size_type j_exp = 0, size_type h_exp = 0)
        :j_exp(j_exp), h_exp(h_exp){
        }
        
        explicit monomial(Vli const& vli, size_type j_exp = 0, size_type h_exp = 0)// for me !
        :j_exp(j_exp), h_exp(h_exp){
            coeff = vli;
        }
        
        monomial& operator *= (Vli const& c){
            (*this).coeff *= c;
            return (*this);
        }
        
        monomial& operator *= (Vli &c){ // I add the &, do you forget it ?
            (*this).coeff *= c;
            return (*this);            
        }
        
        bool operator=(monomial const& m){
            return (*this).coeff == m.coeff;
        }        
        
        typename Vli::value_type const& operator[](size_type i)const{ // for a serial acces of the VLI element
            return coeff[i];
        }

        typename Vli::value_type & operator[](size_type i){
            return coeff[i];
        }

        size_type j_exp;
        size_type h_exp;
        Vli coeff;                
    };
    
    template<class Vli> // need for the boost equal test
    std::ostream& operator<<(std::ostream& os, monomial<Vli> const& m){
        os << m.coeff;
        return os;
    }
}
#endif //VLI_MONOME_H
