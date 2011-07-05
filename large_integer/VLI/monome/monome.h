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
        typedef typename Vli::size_type size_type;
         
        /**
         * Constructor: Creates a monomial 1*J^j_exp*h^h_exp
         */
        explicit monomial(size_type j_exp = 0, size_type h_exp = 0)
        :j_exp(j_exp), h_exp(h_exp){
        }
        
        explicit monomial(Vli const& vli, size_type j_exp = 0, size_type h_exp = 0)// for me !
        :j_exp(j_exp), h_exp(h_exp), coeff(vli){
        }
        
        monomial& operator *= (Vli const& c){
            (*this).coeff *= c;
            return (*this);
        }
        
        monomial& operator *= (int c){
            (*this).coeff *= c;
            return (*this);
        }

        void print(std::ostream& os) const
        {
//            if(coeff > 0)
//                os<<"+";
            os<<coeff<<"*J^"<<j_exp<<"*h^"<<h_exp;
        }

        bool operator == (monomial const& m) const{
            return ((*this).j_exp == m.j_exp) && ((*this).h_exp == m.h_exp) && ((*this).coeff == m.coeff);
        }
        
        typename Vli::value_type const& operator[](size_type i)const{ // for a serial acces of the VLI element
            return coeff[i];
        }

        typename Vli::value_type & operator[](size_type i){
            return coeff[i];
        }
        
        typename Vli::value_type*  p(){
            return coeff.p();
        }

        typename Vli::value_type const* p() const{
            return coeff.p();
        }
        
        size_type j_exp;
        size_type h_exp;
        Vli coeff;
    };
    
    template<class Vli>
    std::ostream& operator<<(std::ostream& os, monomial<Vli> const& m){
        m.print(os);
        return os;
    }
}
#endif //VLI_MONOME_H
