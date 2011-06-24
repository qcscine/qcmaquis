/*
 *  monome.h
 *  vli
 *
 *  Created by Tim Ewart (timothee.ewart@unige.ch) and  Andreas Hehn (hehn@phys.ethz.ch) on 18.03.11.
 *  Copyright 2011 University of Geneva and Eidgenössische Technische Hochschule Züric. All rights reserved.
 *
 */

/*
 * due to the structure of VLI and VECTOR_VLI, be carreful on the operator = 
 */


#ifndef VLI_MONOME_H
#define VLI_MONOME_H

#include <ostream>
#include <cmath>


#define POLYNOMIAL_MAX_ORDER 4 // for testing

namespace vli
{	
	template <class VLI>
    struct monomial
    {
        typedef typename VLI::size_type size_type;	        
        typedef typename VLI::value_type value_type;	
        size_type j_exp;
        size_type h_exp;
        VLI* coeff;
        
        /**
         * Constructor: Creates a monomial 1*J^j_exp*h^h_exp
         */
        explicit monomial(size_type j_exp = 0, size_type h_exp = 0):j_exp(j_exp), h_exp(h_exp){
            coeff = new VLI();
        }
        
        ~monomial(){
            delete coeff;
        }
        
        void destroy(){
            this->~monomial(); //in case of
        }

        monomial& operator *= (VLI const& c);
        monomial& operator *= (value_type c);
        
        
    };
	
	template<class VLI_VECTOR>
	class polynomial{
	public:
		typedef typename VLI_VECTOR::size_type size_type;	
		typedef typename VLI_VECTOR::BaseInt BaseInt;	
        enum {max_order = POLYNOMIAL_MAX_ORDER};
        
        polynomial(){
            coeffs = new VLI_VECTOR(max_order*max_order); // by default everything is set to 0
        }

        polynomial(const polynomial& p){
            coeffs = new VLI_VECTOR(p.coeffs);
        }

        ~polynomial(){
            delete coeffs;
        }
 
        void destroy(){
            this->~polynomial();  //in case of
        }
        
        template<class VLI>
        friend polynomial operator * (const polynomial & p, const monomial<VLI> & m);
  
        polynomial& operator = (polynomial p);

        
        void copy_from_monome(BaseInt* p,const monomial<BaseInt>& m);
        
        friend void swap(polynomial& p1, polynomial& p2);
        inline BaseInt* operator ()(size_type j_exp, size_type h_exp);
        inline const BaseInt* operator ()(size_type j_exp, size_type h_exp) const ; //need ?
        
    private:
        VLI_VECTOR* coeffs;
    };
}

#include "monome/monome.hpp"

#endif //VLI_MONOME_H
