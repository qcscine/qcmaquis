/*
 *  vli.h
 *  vli_cpu_gpu
 *
 *  Created by Tim Ewart on 18.03.11.
 *  Copyright 2011 University of Geneva. All rights reserved.
 *
 */


#ifndef __VLI__
#define __VLI__

#include "definition.h"

namespace vli
{


	template<class T>
	class vli_cpu : public bvli
	{
	public:
		
		vli_cpu()
		{
			size_ = SIZE_BITS/(8*sizeof(T)); 
			data_.resize(size_);
		}	
		
		vli_cpu(T num)
		{
			size_ = SIZE_BITS/(8*sizeof(T)); 
			data_.resize(size_);
			data_[0] = num;
		}		
		
		~vli_cpu()
		{
			data_.erase(data_.begin(),data_.end());
		}
		
		size_int size() const
		{
			return size_;
		}
		
		T &operator[](T i) 
		{
			return data_[i];
		}

		T const & operator[](T i) const 
		{
			return data_[i];
		}
		
		typename std::vector<T>::const_iterator begin()
		{
			typename std::vector<T>::const_iterator it= data_.begin();
			return it;
		}

		typename std::vector<T>::const_iterator end()
		{
			typename std::vector<T>::const_iterator it= data_.end();
			return it;
		}
		
			
	private:
		std::vector<T> data_;
		size_int size_;
		
	};

	template<typename T>
	std::ostream& operator<< (std::ostream& os,  vli_cpu<T> & vli)
	{
		for(typename std::vector<T>::const_iterator it= vli.end(); it != vli.begin(); it--)
			std::cout << *it << " " ;
		os<<std::endl; 
		
		return os;
	}
	



}

#endif
