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
 
CSMatrix& CSMatrix::operator=(const CSMatrix& SMatrix)
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
			_pSData = new CSData(m,n);
			_pSData->PlusRef();	
		}
		
		_pSData->MinusRef();
		_pSData=SMatrix._pSData;
		_pSData->PlusRef();
	}
	
	return *this;
}

double& CSMatrix::operator()(const INTEGER i,const INTEGER j)
{
	if (*(_pSData->_pRefCount) > 1)
	{
		CSData* pSData = new CSData(_pSData);
		_pSData->MinusRef();
		_pSData = pSData;
		_pSData->PlusRef();
	}
		
	return _pSData->_pData->_pArray[i*_m+j];
}


/* ----------------------------------------------------------- Vectors */



CSVector& CSVector::operator=(const CSVector& SVector)
{
	if (_pSData != SVector._pSData)
	{
		if (GetnNumSize() == 0)
		{
			INTEGER mn = SVector.GetnNumSize();
			SetnNumSize(mn);
			_pSData = new CSData(mn);
			_pSData->PlusRef();	
		}
		
		_pSData->MinusRef();
		_pSData=SVector._pSData;
		_pSData->PlusRef();
	}
	
	return *this;
}

double& CSVector::operator()(const INTEGER i)
{
	if (*(_pSData->_pRefCount) > 1)
	{
		CSData* pSData = new CSData(_pSData);
		_pSData->MinusRef();
		_pSData = pSData;
		_pSData->PlusRef();
	}
		
	return _pSData->_pData->_pArray[i];
}

