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
*/

/**
* size of the vli
*/
#define SIZE_VLI 4
/**
* size of the polynomial
*/
#define SIZE_POLY_VLI 4
/**
* size of the vector
*/
#define SIZE_VECTOR_VLI 4
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

    struct size_vector_vli    
    {
        enum { value = SIZE_VECTOR_VLI};
    };

    } // end namespace detail
} // end namespace vli
#endif  //VLI_SIZE_PARAM_HPP


