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
//#define LOG_BASE			 0x1E // 30
//#define BASE_MINUS1			 0x3FFFFFFF
#define LOG_BASE			 2 // 30
#define BASE				 4
#define BASE_MINUS2			 2
#define MINUS_BASE_PLUS2	 -2


typedef std::size_t size_type;

template <class T>
void addition(T A, T B, T C,int num_integer, int ld);  

#endif