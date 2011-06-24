/*
 *  monome.h
 *  vli
 *
 *  Created by Tim Ewart (timothee.ewart@unige.ch) and  Andreas Hehn (hehn@phys.ethz.ch) on 24.06.11.
 *  Copyright 2011 University of Geneva and Eidgenössische Technische Hochschule ZüricZ. All rights reserve.
 *
 */



/**
 * Multiplication of a polynomial with a monomial
 */

#ifndef VLI_MONOME_HPP
#define VLI_MONOME_HPP

namespace vli
{	
/*************************** MONOME *************************************/
    
template<class VLI>
monomial<VLI>& monomial<VLI>::operator *= (VLI const& c)
{
    coeff *= c;
    return *this;
}

template<class VLI>    
monomial<VLI>& monomial<VLI>::operator *= (typename VLI::value_type c)    
{
    assert(false); // TO DO !!
    coeff *= c;
    return *this;
}
    

/*************************** POLYNOME ***********************************/
    
template <class BaseInt> //CPU specialization 
polynomial<vli_vector_cpu<BaseInt> > operator* (polynomial<vli_vector_cpu<BaseInt> > const& p, monomial<vli_cpu<BaseInt> > const& m)
{
    typedef typename monomial<vli_cpu<BaseInt> >::size_type size_type;
    typedef typename monomial<vli_cpu<BaseInt> >::value_type value_type;
        
    polynomial<vli_vector_cpu<BaseInt> > r;
    vli_cpu<BaseInt> c;
    size_type max_order = r.max_order;
    for(std::size_t je = 0; je < max_order-m.j_exp; ++je){
        for(std::size_t he = 0; he < max_order-m.h_exp; ++he){
            /******  r(je+m.j_exp,he+m.h_exp) = p(je,he) * m.coeff; ORIGINAL ******/
            memcpy((void*)&c[0], (void*)p(je,he),full_size_single_vli);
            multiplies_assign(c,m.coeff);
            memcpy((void*)p(je,he),(void*)&c[0],full_size_single_vli);
        }
    }
    return r;
}
    
/**
* Assignment operator
*/
template <class VLI_VECTOR>
polynomial<VLI_VECTOR>& polynomial<VLI_VECTOR>::operator = (polynomial p)
{
    swap(*this,p);
    return *this;
}
 /*
template <class VLI_VECTOR>
void  polynomial<VLI_VECTOR>::copy_from_monome(typename VLI_VECTOR::BaseInt* p,const  monomial<vli_cpu<VLI_VECTOR::BaseInt> >& m)
{
    memcpy((void*)this->p, (void*)&monomial.coeff[0], full_size_single_vli);
}
   */ 
/**
* Swap function
*/
template <class VLI_VECTOR>
void swap(polynomial<VLI_VECTOR>& p1, polynomial<VLI_VECTOR>& p2)    
{
    using std::swap;
    swap(p1.coeffs,p2.coeffs);
}    
    
/** 
* Access coefficient of monomial J^j_exp*h^h_exp, here we have a pointer on the first element of VLI
*/
template<class VLI_VECTOR>
inline typename VLI_VECTOR::BaseInt* polynomial<VLI_VECTOR>::operator ()(typename VLI_VECTOR::size_type j_exp, typename VLI_VECTOR::size_type h_exp) 
{
    assert(j_exp < max_order);
    assert(h_exp < max_order);
    return coeffs[j_exp*max_order+h_exp];
}
    
template<class VLI_VECTOR>
inline const typename VLI_VECTOR::BaseInt* polynomial<VLI_VECTOR>::operator ()(typename VLI_VECTOR::size_type j_exp, typename VLI_VECTOR::size_type h_exp) const 
{
    assert(j_exp < max_order);
    assert(h_exp < max_order);
    return coeffs[j_exp*max_order+h_exp];
}
    

}//end namespace

#endif
