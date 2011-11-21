/*
 *  list_functions.h
 *  dmrg
 *
 *  Created by Tim Ewart on 29.01.11.
 *  Copyright 2011 Université de Genève. All rights reserved.
 *
 */

/**
	the two template arguments are the p_ pointer inside gpu_matrix

 */
/*
template <class T>
void transpose(T* A, T* B, std::size_t num_rows, std::size_t num_columns, std::size_t ld);  
*/


template <class T>
void transpose(T A, T B, std::size_t   num_rows, std::size_t   numcolumns, std::size_t   ld);  

template <class T>
void swap_rows(T A, std::size_t num_rows, std::size_t num_columns, std::size_t ld , std::size_t i1, std::size_t i2);

template <class T>
void swap_columns(T A, std::size_t num_rows, std::size_t num_columns, std::size_t ld , std::size_t i1, std::size_t i2);



//void swap_columns(size_type j1, size_type j2);