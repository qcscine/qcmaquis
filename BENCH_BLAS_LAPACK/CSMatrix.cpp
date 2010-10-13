/*
 *  CSMatrix.cpp
 *  SMART_MATRIX
 *
 *  Created by Tim Ewart on 29.09.10.
 *  Copyright 2010 University of Geneva. All rights reserved.
 *
 */
#include <iostream>

#ifdef VECTOR_STL 
	#include <vector>
#endif

#include "CSMatrix.h"

/* ----------------------------------------------------------- Matrices */
template <class Datatype>
CSMatrix<Datatype>& CSMatrix<Datatype>::operator=(const CSMatrix <Datatype>& SMatrix)
{
	if (_pSData != SMatrix._pSData)
	{
		if (GetnNumSize() == 0)
		{
			INTEGER m = SMatrix.GetnNumRow();
			INTEGER n = SMatrix.GetnNumCol();
			SetnNumRow(m);
			SetnNumCol(n);
			SetnNumSize(m*n);
			_pSData = new CSData<Datatype>(m,n);
			_pSData->PlusRef();	
		}
		
		_pSData->MinusRef();
		_pSData=SMatrix._pSData;
		_pSData->PlusRef();
	}
	
	return *this;
}

template <class Datatype>
Datatype& CSMatrix<Datatype>::operator()(const INTEGER i,const INTEGER j)
{
	if (*(_pSData->_pRefCount) > 1)
	{
		CSData<Datatype>* pSData = new CSData<Datatype>(_pSData);
		_pSData->MinusRef();
		_pSData = pSData;
		_pSData->PlusRef();
	}
		
	return _pSData->_pData->_pArray[i*_m+j];
}


/* ----------------------------------------------------------- Vectors */


template <class Datatype>
CSVector<Datatype>& CSVector<Datatype>::operator=(const CSVector<Datatype>& SVector)
{
	if (_pSData != SVector._pSData)
	{
		if (GetnNumSize() == 0)
		{
			INTEGER mn = SVector.GetnNumSize();
			SetnNumSize(mn);
			_pSData = new CSData<Datatype>(mn);
			_pSData->PlusRef();	
		}
		
		_pSData->MinusRef();
		_pSData=SVector._pSData;
		_pSData->PlusRef();
	}
	
	return *this;
}

template <class Datatype>
Datatype& CSVector<Datatype>::operator()(const INTEGER i)
{
	if (*(_pSData->_pRefCount) > 1)
	{
		CSData<Datatype>* pSData = new CSData<Datatype>(_pSData);
		_pSData->MinusRef();
		_pSData = pSData;
		_pSData->PlusRef();
	}
		
	return _pSData->_pData->_pArray[i];
}

