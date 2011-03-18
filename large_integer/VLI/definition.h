/*
 *  definition.h
 *  Untitled
 *
 *  Created by Tim Ewart on 21.02.11.
 *  Copyright 2011 University of Geneva. All rights reserved.
 *
 */

#ifndef __DEFINITION__
#define __DEFINITION__

#define NUM 16;

// BASE 255
/*
#define LOG_BASE			 0x8   // 8 
#define LOG_BASE_HALF		 0x4   // 4
#define BASE				 0x100 // 256
#define BASE_HALF			 0x10  // 16
#define BASE_MINUS2			 0xFE  // 254   
#define MINUS_BASE_PLUS2	 -254
#define MASK_DOWN			 0xF
#define MASK_UP				 0xF0
*/

//BASE 2**30
#define LOG_BASE			 30   
#define LOG_BASE_HALF		 15   
#define BASE				 0x7FFFFFFF 
#define BASE_HALF			 0x3FFFFFFF 
#define BASE_MINUS2			 0x7FFFFFFD
#define MINUS_BASE_PLUS2	 -0x7FFFFFFD
#define MASK_DOWN			 0xFFFF
#define MASK_UP				 0xFFFF0000


typedef std::size_t size_int;
typedef std::size_t size_type;

class bvli
{
public:
	bvli(){};
	virtual ~bvli(){};	
	
};


template <class T>
void addition(T A, T B, T C,int num_integer, int ld);  

template <class T>
void multiplication(T A, T B, T C,int num_integer, int ld);  





#endif