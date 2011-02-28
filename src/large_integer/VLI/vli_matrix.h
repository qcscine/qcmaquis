/*
 *  sparse_matrix.h
 *  Untitled
 *
 *  Created by Tim Ewart on 21.02.11.
 *  Copyright 2011 University of Geneva. All rights reserved.
 *
 */

#ifndef __SPARSE_MATRIX__
#define __SPARSE_MATRIX__

#include <vector>
#include <iterator> 
#include <stdlib.h>
#include <assert.h>
#include <iostream>

#include "definition.h"


/**
vli means very long integer 
*/

namespace vli
{

template <typename T>
class vli_matrix
{
public:
	
/**
 length_vli = 4 -> 120 bit integer 
*/	
	vli_matrix(size_type numcols, size_type length_vli = 4 )
	{
		numrows_ = length_vli;
		numcols_ = numcols;
		ld_ = length_vli;
		length_ = numcols*length_vli;
		data_.resize(length_);
	}
	
	
	vli_matrix(vli_matrix &  matrix)
	{
		numrows_ = matrix.num_rows();
		numcols_ = matrix.num_columns();
		ld_ = matrix.stride2();
		length_ = numcols_*numrows_;
		data_.resize(length_);
	}
	
	/**
	To have the compatibility with my gpu matrix class
	I need the followings five functions 
	*/
	inline T&  operator()( size_type i,  size_type j) 
	{		
		return data_[i + numrows_*j];
	}
	
	inline T const&  operator()( size_type i,  size_type j) const
	{		
		return data_[i + numrows_*j];
	}
	
	inline const size_type num_columns() const
	{
		return numcols_;
	}

	inline const size_type  num_rows() const 	
	{
		return numrows_;
	}
	
	inline const size_type  stride2() const	
	{
		return ld_;
	}
	
	
	void init_random()
	{
		typename std::vector<T>::iterator it;
		for(it = data_.begin(); it != data_.end() ; it++ )
			*it = rand()%16;		
	
	}
		
private:

	size_type numcols_;
	size_type numrows_;
	size_type ld_;
	std::vector<T> data_ ;
	size_type length_;
};

	
template<typename T>
std::ostream& operator<< (std::ostream& os,  vli_matrix<T> & matrix_cpu)
{
		for (size_type i = 0 ; i < matrix_cpu.num_rows(); i++)
		{
			for(size_type j = 0 ; j < matrix_cpu.num_columns() ; j++)
			{
				os << matrix_cpu(i,j) << " " ;
			}
			os << std::endl;
		}

	return os;
}

template <typename T>
void addition_Avizienis_kernel_cpu(vli::vli_matrix<T> const & x, vli::vli_matrix<T> const &  y , vli::vli_matrix<T>  &  z);

	
	
}
#endif


