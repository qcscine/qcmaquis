
typedef double __attribute__((aligned(16))) SSE_PTRDOUBLE; //Alignement si GPU un jour (version gcc)

class CMatrix
{
public:
	
	CMatrix( int nRows,  int nColumns):_nRows(nRows),_nColumns(nColumns)
	{
		_ppM = new SSE_PTRDOUBLE*[nColumns];
		_Contiguous_Array = new double[nRows*nColumns]; 
		int j=0;
		for( int i=0;i<nColumns;i++)
		{
			_ppM[i] = &_Contiguous_Array[0]+j;
			j=j+nRows;
		}
	};
	
	~CMatrix()
	{		
	 
		delete [] _Contiguous_Array;
		delete [] _ppM;
	}
	
	// Attention alignement en ligne RowMajor, or lapack veut des columnMajor, astuce à la Tim
	/*
	 CBLASS peut accepter des alignements lignes (C) et colonnes (F77) grace aux variables d'ajustement  CblasColMajor et CblasRowMajor
	 CLAPACK n'accepte que des alignements colonnes 
	 conclusion astuce pour le remplissage du vecteur de la matrice
	 */
	double& operator()( int i, int j) const // référence constante on peut l'utiliser en left et right value
	{
		return _ppM[j][i]; // inversion des indices i et j pour avoir un alignement en colonnes
	}
	
	
	CMatrix& operator=(const CMatrix& A)
	{
		for(int i=0;i<_nColumns;i++)
		{
			for(int j=0;j<_nRows;j++)
			{
				_ppM[i][j] = A._ppM[i][j];
			}
		}	
		return *this;
	}
	
	
	
	void Print()
	{
		for(int i=0;i<_nRows;i++)
		{
			for(int j=0;j<_nColumns;j++)
			{
				std::cout << _ppM[j][i] << " " ;
			}
			std::cout << std::endl;
		}
	}
	
	int  GetnRows() const
	{
		return _nRows;
	}
	
	int GetnCols() const
	{
		return _nColumns;
	}
	
	void InitRandom();
	
private:
	double** _ppM;
	double* _Contiguous_Array;
	int  _nRows;
	int  _nColumns;
};



