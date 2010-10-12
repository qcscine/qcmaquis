/*
 *  CSMatrix.h
 *  SMART_MATRIX
 *
 *  Created by Tim Ewart on 29.09.10.
 *  Copyright 2010 University of Geneva. All rights reserved.
 *
 */


typedef double __attribute__((aligned(16))) SSE_PTRDOUBLE; /*GCC Only, QuadWord alignement i.e. 16 byte for SIMD and GPU*/
typedef long int INTEGER; /*Avoid fortran(int 8 bytes) and C(int 4 bytes) conflicts*/

// Is it just a data container, it is a basic STL vector !
class CArray
{
public:
	CArray(INTEGER n):_nSize(n)
	{
#ifdef VECTOR_STL
		_pArray.resize(n);
#else
		_pArray = new double[n];
#endif
	}
	
	~CArray()
	{
#ifdef VECTOR_STL
		_pArray.clear();
#else		
		delete [] _pArray;
#endif			
	}
	
	INTEGER Size() const {return _nSize;} 

	
#ifdef VECTOR_STL
	std::vector <double> _pArray;
#else
	SSE_PTRDOUBLE * _pArray;
#endif	
	INTEGER _nSize;
};

//Smart Pointer class
class CSData
{
	public:	
	//For Matrices
	CSData(INTEGER m, INTEGER n)
	{
		_pData = new CArray(m*n); 
		_pRefCount = new INTEGER(0);
	}
	
	//For Vectors
	CSData(INTEGER mn)
	{
		_pData = new CArray(mn); 
		_pRefCount = new INTEGER(0);
	}
	
	//For Vectors and Matrices
	CSData(CSData* pSData)
	{
		INTEGER Size = pSData->_pData->Size();
		_pData = new CArray(Size);
		
		for(INTEGER i = 0; i< Size; i++)
			_pData->_pArray[i] = pSData->_pData->_pArray[i];
		
		_pRefCount = new INTEGER(0);
		
		
	};
	
	~CSData()
	{
		std::cout << "delete CSData" << std::endl;
		if (_pRefCount != 0)
		{
			(*_pRefCount)--;
			
			if ((*_pRefCount) <= 0)
			{
				delete _pRefCount;
				delete _pData;
			}
		}

	};
	
	void PlusRef()
	{
		(*_pRefCount)++;
	}
	
	void MinusRef()
	{
		if ((*_pRefCount)-- == 0)
		{
			delete this;
		}
	}

	INTEGER* _pRefCount; //counter of reference	
	CArray * _pData; //data	
	
};

//Matrix class
class CSMatrix
{
public:
	
	CSMatrix(INTEGER m = 0, INTEGER n = 0):_m(m),_n(n),_mn(m*n)
	{
		if (_mn != 0) 
		{
			//smart pointer for data
			_pSData = new CSData(m,n);
			_pSData->PlusRef();		
		}
	}
	
	
	~CSMatrix()
	{
		std::cout << "delete matrice :"<< std::endl; //just to be sure, to remove
		delete _pSData;
	}
		
	void Init(INTEGER a)
	{
		
		srand(2);
		for(INTEGER i=0;i < _mn; i++)
		{
			_pSData->_pData->_pArray[i] = rand()*3.14;
		}
		
	}
	
	void Print()
	{
		for(int i=0;i < GetnNumRow()*GetnNumCol(); i++)
			std::cout <<  _pSData->_pData->_pArray[i] << std::endl;
	
	}
	
	INTEGER GetnNumRow() const
	{
		return _m; 
	}
	
	INTEGER GetnNumCol() const
	{
		return _n; 
	}
	
	INTEGER GetnNumSize() const
	{
		return _mn;
	}
	
	void SetnNumRow(INTEGER m) 
	{
		 _m = m; 
	}
	
	void SetnNumCol(INTEGER n) 
	{
		 _n = n; 
	}
	
	void SetnNumSize(INTEGER mn) 
	{
		 _mn = mn;
	}
	
	
	
	
	CSMatrix& operator=(const CSMatrix& SMatrix);
	double& operator()(const INTEGER i, const INTEGER j);
	
	
private:

	CSData *_pSData;
	INTEGER _m, _n, _mn;
};

//Vector class
class CSVector
{
public:
	
	CSVector(INTEGER mn):_mn(mn)
	{
		if (_mn != 0) 
		{
			_pSData = new CSData(mn);
			_pSData->PlusRef();
		}
	}
	
	
	~CSVector()
	{
		std::cout << "delete vector :"<< std::endl; //just to be sure, to remove
		delete _pSData;
	}
	
	void Init(INTEGER a)
	{
		for(INTEGER i=0;i < _mn; i++)
		{
			_pSData->_pData->_pArray[i] = static_cast<double>(a+i);
		}
		
	}
	
	void Print()
	{
		for(int i=0;i < _mn ;i++)
			std::cout <<  _pSData->_pData->_pArray[i] << std::endl;
		
	}
		
	INTEGER GetnNumSize() const
	{
		return _mn;
	}
	
	void SetnNumSize(INTEGER mn) 
	{
		 _mn = mn;
	}
	
	CSVector& operator=(const CSVector& SMatrix);
	double& operator()(const INTEGER i);
	
	
private:
	
	CSData *_pSData;
	INTEGER _mn;
	
	
};

