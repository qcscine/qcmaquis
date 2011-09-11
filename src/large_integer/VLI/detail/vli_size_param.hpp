//
//  vli_size_param.h
//  vli
//
//  Created by Timothée Ewart on 23/08/11.
//  Copyright 2011 IBM. All rights reserved.
//

/**
*
*
*/

#ifndef VLI_SIZE_PARAM_HPP
#define VLI_SIZE_PARAM_HPP

/** BEGINNING FIX YOUR PARAMETERS
* Fix your parameters before compiling !!!!!!!
* Parameters are fixed into cmake
*/

/**
* size of the polynomial
*/
#define SIZE_POLY_VLI POLYNOMIAL_MAX_ORDER
/**
* size of the vector
*/
#define SIZE_VECTOR_VLI VECTOR_MAX_ORDER
/**
* type of the vli
*/
typedef unsigned long int type;

/** 
*END FIX YOUR PARAMETERS 
*/

/**
*   encapsuled parameters
*/
namespace vli
{
    namespace detail
    {

    struct size_poly_vli
    {
        enum { value = SIZE_POLY_VLI};
    };

    } // end namespace detail
} // end namespace vli
#endif  //VLI_SIZE_PARAM_HPP


