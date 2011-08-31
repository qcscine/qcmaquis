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
* size of the vli
*/
#define SIZE_VLI VLI_MAX_ORDER
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
typedef unsigned int type;

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

    struct type_vli
    {
        typedef type BaseInt ;
    };

    struct size_vli
    { 
        enum { value = SIZE_VLI};
    };

    struct size_poly_vli
    {
        enum { value = SIZE_POLY_VLI};
    };

    // Presently not used, because the vectors are dynamics
    struct size_vector_vli    
    {
        enum { value = SIZE_VECTOR_VLI};
    };

    } // end namespace detail
} // end namespace vli
#endif  //VLI_SIZE_PARAM_HPP


