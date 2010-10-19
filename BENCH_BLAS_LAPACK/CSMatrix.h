/*
 *  CSMatrix.h
 *  SMART_MATRIX
 *
 *  Created by Tim Ewart on 29.09.10.
 *  Copyright 2010 University of Geneva. All rights reserved.
 *
 */
typedef double __attribute__((aligned(16))) SSE_PTRDOUBLE; /*GCC Only, QuadWord alignement i.e. 16 byte for SIMD and GPU*/
typedef  int INTEGER; /*Avoid fortran(int 8 bytes) and C(int 4 bytes) conflicts*/

// Is it just a data container, it is a basic STL vector !

template <class DataType>
class CArray
{
public:
	CArray(INTEGER n):_nSize(n)
	{
#ifdef VECTOR_STL
		_pArray.resize(n);
#else
		_pArray = new DataType[n];
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
	std::vector <DataType> _pArray;
#else
//	SSE_PTRDOUBLE * _pArray; 
	DataType * _pArray; //Allign on qudWord
#endif	
	INTEGER _nSize;
};

//Smart Pointer class
template <class DataType>
class CSData
{
	public:	
	//For Matrices
	CSData(INTEGER m, INTEGER n)
	{
		_pData = new CArray<DataType> (m*n); 
		_pRefCount = new INTEGER(0);
	}
	
	//For Vectors
	CSData(INTEGER mn)
	{
		_pData = new CArray<DataType>(mn); 
		_pRefCount = new INTEGER(0);
	}
	
	//For Vectors and Matrices
	CSData(CSData* pSData)
	{
		INTEGER Size = pSData->_pData->Size();
		_pData = new CArray<DataType>(Size);
		
		for(INTEGER i = 0; i< Size; i++)
			_pData->_pArray[i] = pSData->_pData->_pArray[i];
		
		_pRefCount = new INTEGER(0);
		
		
	};
	
	~CSData()
	{
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
	CArray<DataType> * _pData; //data	
	
};

//Matrix class
template <class DataType> 
class CSMatrix
{
public:
	
	CSMatrix(INTEGER m = 0, INTEGER n = 0):_m(m),_n(n),_mn(m*n)
	{
		if (_mn != 0) 
		{
			//smart pointer for data
			_pSData = new CSData<DataType>(m,n);
			srand(2);
			for(INTEGER i=0;i < _mn; i++)
			{
				
#ifdef 	DOUBLE			
				_pSData->_pData->_pArray[i] = rand()*3.14;
#endif
#ifdef COMPLEX
				_pSData->_pData->_pArray[i].real() = rand()/3.14;
				_pSData->_pData->_pArray[i].imag() = rand()*3.14;

				
#endif				
				
			}
			_pSData->PlusRef();		
		}
	}
	
	
	~CSMatrix()
	{
		delete _pSData;
	}
		

	
	void Print()
	{
		
		for(int i =0 ; i< _n;i++)
		{
			for(int j=0 ; j< _m;j++)
			{
				std::cout <<  _pSData->_pData->_pArray[i*_m+j] << " " ;
			}	
			std::cout << std::endl;
		}
		
		
	//	for(int i=0;i < GetnNumRow()*GetnNumCol(); i++)
	//		std::cout <<  _pSData->_pData->_pArray[i] << std::endl;
	
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
	
	
	
	
	CSMatrix<DataType>& operator=(const CSMatrix<DataType>& SMatrix);
	DataType& operator()(const INTEGER i, const INTEGER j);
	
	
private:

	CSData<DataType> *_pSData;
	INTEGER _m, _n, _mn;
};

//Vector class
template <class DataType> 
class CSVector
{
public:
	
	CSVector(INTEGER mn):_mn(mn)
	{
		if (_mn != 0) 
		{
			_pSData = new CSData<DataType>(mn);
			_pSData->PlusRef();
		}
	}
	
	
	~CSVector()
	{
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
	DataType& operator()(const INTEGER i);
	
	
private:
	
	CSData<DataType> *_pSData;
	INTEGER _mn;	
	
};

