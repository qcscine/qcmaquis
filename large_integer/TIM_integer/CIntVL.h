/*
 *  CVlint.h
 *  512BITS
 *
 *  Created by Tim Ewart on 02.11.10.
 *  Copyright 2010 University of Geneva. All rights reserved.
 *
 */

template <class Type_Container, class Type>
class Clint
{
public:
	
	Clint():m_nNBit(128){}
	
	Clint(const Clint<Type_Container, Type>& lintA):m_nNBit(128),m_Data(lintA.m_Data){}
	
	~Clint(){}
	
	
	size_t GetSizeBits() const
	{
		return m_nNBit;
	}	
	
	long long int &operator[](size_t nIndex)
	{
		return m_Data[nIndex];
	}
	
	Clint<Type_Container, Type>& operator=(Type A)
	{
		m_Data[0] = static_cast<unsigned long long int> (A);
		return *this;
	}
	
	
	Clint<Type_Container, Type>& operator=(const Clint<Type_Container,Type> & lintA)
	{
		m_Data = lintA.m_Data;
		return *this;
	}
	
	
	friend Clint<Type_Container, Type> operator+(const Clint<Type_Container , Type> & lintA, const Clint<Type_Container, Type> & lintB)
	{
		Clint<Type_Container,Type> lintC;
		
		lintC.m_Data = lintA.m_Data + lintB.m_Data;
		
		return lintC;
	}
		
	
private:
	
	Type_Container m_Data;
	std::size_t m_nNBit; //Presently no interest
	
};
