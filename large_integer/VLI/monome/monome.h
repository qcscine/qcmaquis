/*
 *  monome.h
 *  monome
 *
 *  Created by Tim Ewart on 18.03.11.
 *  Copyright 2011 University of Geneva. All rights reserved.
 *
 */


#ifndef VLI_MONOME_H
#define VLI_MONOME_HP

#include <ostream>
#include <cmath>
#include <vector>

namespace vli
{	
	template<class VLI>
	class vli_site_monome{
	public:
        typedef typename VLI::size_type size_type;		

		explicit vli_site_monome(size_type exponent = 0):exponent_(exponent){
			mono_ = new VLI();
	    }
				
		vli_site_monome(const vli_site_monome & monome_site){
    		this->mono_ = new VLI(*monome_site.mono_);
			this->exponent_ = monome_site.exponent_;
		}
				
		~vli_site_monome(){
			delete mono_;
		}
		
		vli_site_monome & operator=(vli_site_monome site_monome){
			this->swap(site_monome);
			return *this;
		}

		void swap(vli_site_monome & site_monome){
			std::swap(this->mono_, site_monome.mono_);
			std::swap(this->exponent_, site_monome.exponent_);
		} 
		
		size_type & exponent(){
		   	return exponent_;
		}
		
	private:
		VLI* mono_;
	    size_type exponent_;
	};
	
	template<class VLI>
	class vli_monome{
	public:
		typedef typename VLI::size_type size_type;	
	
		explicit vli_monome(){
		    vec_.push_back(new vli_site_monome<VLI>());
		}
				
		vli_monome(vli_monome & monome){ // no const ...
			for(size_type i(0);i<monome.size();i++){
			    this->vec_.push_back(new vli_site_monome<VLI>(*monome(i)));
			}
		}
		
		vli_monome & operator=(vli_monome monome){
            swap(monome);
			return *this;
		}

		void swap(vli_monome & monome){
			std::swap(this->vec_, monome.vec_);	
		}
				
		template<class T> // for iterator
		vli_site_monome<VLI>* operator()(T t){
			return vec_[t];
		}
		
		size_type & operator[](size_type i){
		    return vec_[0]->exponent();	
		}
		
		const size_type size(){
			return vec_.size();
		}
		
	private:
		std::vector<vli_site_monome<VLI>* > vec_;
	};
}

#endif //VLI_MONOME_H
