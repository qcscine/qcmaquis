/*
 *  monome.h
 *  monome
 *
 *  Created by Tim Ewart on 18.03.11.
 *  Copyright 2011 University of Geneva. All rights reserved.
 *
 */


#ifndef VLI_MONOME_H
#define VLI_MONOME_H

#include <ostream>
#include <cmath>
#include <vector>

#define MAX_EXPONENT 0x2 //today ...

namespace vli
{	
	template <class VLI>
    struct monomial_gpu
    {
        typedef typename VLI::size_type size_type;	
        size_type j_exp;
        size_type h_exp;
        VLI coeff;
        
        /**
         * Constructor: Creates a monomial 1*J^j_exp*h^h_exp
         */
        explicit monomial(size_type j_exp = 0, size_type h_exp = 0)
        : j_exp(j_exp), h_exp(h_exp), coeff(1)
        {
        }
        
        monomial& operator *= (CoeffType const& c)
        {
            coeff *= c;
            return *this;
        }
        
        monomial& operator *= (int c)
        {
            coeff *= c;
            return *this;
        }
    };
	
	template<class VLI>
	class vli_polynomial_gpu{
	public:
		typedef typename VLI::size_type size_type;	
       // enum {max_order = POLYNOMIAL_MAX_ORDER};
        enum {max_order = 10};

        
        vli_polynomial(){
            coeffs_.resize(max_order*max_order);

        }
        
    private:
        std::vector<vli_monome<VLI>* > coeffs_;
    };
}

#endif //VLI_MONOME_H
