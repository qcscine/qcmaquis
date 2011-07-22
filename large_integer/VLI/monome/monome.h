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
        :j_exp_(j_exp), h_exp_(h_exp){
        }
        
        explicit monomial(Vli const& vli, size_type j_exp = 0, size_type h_exp = 0)// for me !
        :j_exp_(j_exp), h_exp_(h_exp), coeff_(vli){
        }
        
        monomial& operator *= (Vli const& c){
            coeff_ *= c;
            return (*this);
        }
        
        monomial& operator *= (int c){
            coeff_ *= c;
            return (*this);
        }

        void print(std::ostream& os) const
        {
            if(coeff_ > 0)
                os<<"+";
            os<<coeff_<<"*J^"<<j_exp_<<"*h^"<<h_exp_;
        }

        bool operator == (monomial const& m) const{
            return (j_exp_ == m.j_exp_) && (h_exp_ == m.h_exp_) && (coeff_ == m.coeff_);
        }
        
        typename Vli::value_type const& operator[](size_type i)const{ // for a serial acces of the VLI element
            return coeff_[i];
        }

        typename Vli::value_type & operator[](size_type i){
            return coeff_[i];
        }
        
        typename Vli::value_type*  p(){
            return coeff_.p();
        }

        typename Vli::value_type const* p() const{
            return coeff_.p();
        }
        
        size_type j_exp_;
        size_type h_exp_;
        Vli coeff_;
    };
    
    template<class Vli>
    std::ostream& operator<<(std::ostream& os, monomial<Vli> const& m){
        m.print(os);
        return os;
    }
    
    template <typename Vli, typename T>
    monomial<Vli> operator * (monomial<Vli> m, T c)
    {
        m*=c;
        return m;
    }

    template <typename Vli, typename T>
    monomial<Vli> operator * (T c, monomial<Vli> m)
    {
        return m*c;
    }

}


#endif //VLI_MONOME_H
