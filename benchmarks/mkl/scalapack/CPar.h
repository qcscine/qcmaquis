
/*
C++ style = reference, pointers are evil, but we used fortran library !!!
*/
extern "C" void blacs_pinfo_( int* mypnum, int* nprocs); 
extern "C" void blacs_gridinfo_( int* nContxt, int* nRows, int* nCols, int* nMyRow, int* nMyCol);
extern "C" void blacs_gridinit_(  int* nContxt, char * order,  int* np_row,  int* np_col);
extern "C" void blacs_gridexit_( int* nContxt);
extern "C" void blacs_exit_( int* error_code);
extern "C" void blacs_get_(int* Contxt, int* what, int* val); 
extern "C" int numroc_(int* nN, int* nB, int* nIProc, int* nISrcProc, int* nNProcs);
extern "C" void pdgesvd_(char *jobu, char *jobvt,
		int *m, int *n,
		double *a, int *ia, int *ja, int *desca,
		double *s,
		double *u, int *iu, int *ju, int *descu,
		double *vt, int *ivt, int *jvt, int *descvt,
		double *work, int *lwork,double *rwork, int *info);


class CGrid
{
public:
	//First Constructor by hand !
	CGrid(int nRowsGrid, int nColsGrid, int nDesiredBlockSiwe):_nRowsGrid(nRowsGrid),_nColsGrid(nColsGrid),
								   _nDesiredBlockSize(nDesiredBlockSiwe),_nMyRow(0),_nMyCol(0){}
	
	//Second Constructor automatic !
	CGrid(int nNumProcs, int nDesiredBlockSiwe):_nMyRow(0),_nMyCol(0)
	{
		_nDesiredBlockSize = nDesiredBlockSiwe;
		//Cast C++ style
		_nRowsGrid = static_cast<int>(floor(sqrt(nNumProcs)));
		
		while (nNumProcs % _nRowsGrid != 0)
			_nRowsGrid--;

		_nColsGrid = nNumProcs/_nRowsGrid;
	};
	

	~CGrid(){}	
	
	void InitBlockSize(int nM, int nN)
	{
		int zero = 0;
	
		if (nM < _nRowsGrid*_nDesiredBlockSize)
		{
			_nBlockSizeRows = nM/_nRowsGrid;
			_nBlockSizeCols = nN/_nColsGrid;
		} 
		else 
		{
			_nBlockSizeRows = _nDesiredBlockSize;
			_nBlockSizeCols = _nDesiredBlockSize;
		}
		
		_nRowDis = numroc_(&nM, &_nBlockSizeRows, &_nMyRow, &zero, &_nRowsGrid);	
		_nColDis = numroc_(&nN, &_nBlockSizeCols, &_nMyCol, &zero, &_nColsGrid);
		
	}
	
	/*
	Bad C style for reference	
	*/	
		
	int* GetnRowsGrid()  {return &_nRowsGrid;}
	int* GetnColsGrid()  {return &_nColsGrid;}	
	int* GetnMyRowsGrid()  {return &_nMyRow;}
	int* GetnMyColsGrid()  {return &_nMyCol;}
	int GetnRowsDisGrid()  {return _nRowDis;}
	int GetnColsDisGrid()  {return _nColDis;}
		
	int GetnBlockSizeRows() const {return _nBlockSizeRows;}
	int GetnBlockSizeCols() const {return _nBlockSizeCols;}
	int& ChangeBlockSize() {return _nDesiredBlockSize;} 
	
private:
	// coordinate Row, Col : specific to a professor
	int _nMyRow, _nMyCol;
	// number of row and col of the grid
	int _nRowsGrid, _nColsGrid;
	//
	int _nDesiredBlockSize;
	//
	int _nBlockSizeRows,_nBlockSizeCols;

	int _nSize;
	// nuber of Rows or Columns of a distributed matrix owned by the process indicated by IPROC.
	int _nRowDis, _nColDis;
};



class CParMPI
{
	public:
	
	CParMPI(int argc, char * argv[])
	{
		MPI_Init(&argc,&argv);
	        MPI_Comm_size(MPI_COMM_WORLD,&_nnumprocsMPI);
       		MPI_Comm_rank(MPI_COMM_WORLD,&_nmyidMPI);		
	}
	
	~CParMPI()
	{
		MPI_Finalize();
	}
	
	
	int  GetMyIdMPI() const {return _nmyidMPI;}
	int  GetNumProcsMPI() const {return _nnumprocsMPI;}

	void Print(CGrid& Grid)
	{
		cout << " num procs " << _nnumprocsMPI << " MyID : " << _nmyidMPI ;
		cout << " , num Rows " << *Grid.GetnRowsGrid() << " num Cols : " << *Grid.GetnColsGrid()  ;		
		cout << " , num myRows " << *Grid.GetnMyRowsGrid() << " num myCols : " << *Grid.GetnMyColsGrid()  << endl;		
	}

	
	private :
	char _param[256];
	int _nmyidMPI,_nnumprocsMPI;
	
};

/*
 CParBLACS derived from CParMPI because we need the mpi initialization first  
*/

class CParBLACS : public CParMPI
{
public:
	CParBLACS(int argc, char * argv []):CParMPI(argc,argv)
	{
		_order[0] = 'R';
		_order[1] = 'o';
		_order[2] = 'w';
		_nContxt = 0;
		_nContinue = 1; //for the destructor	
		_nVal = 0;	
		blacs_pinfo_(&_nmyidBLACS,&_nnumprocsBLACS);	
	 	blacs_get_(&_nContxt,&_nContinue,&_nVal);
	}
	
	~CParBLACS()
	{
		blacs_gridexit_(&_nVal);
		blacs_exit_(&_nContinue);
	}
	
	
	void InitGrid(CGrid& Grid)
	{
		blacs_gridinit_(&_nVal,_order,Grid.GetnRowsGrid(),Grid.GetnColsGrid());
		blacs_gridinfo_(&_nVal,Grid.GetnRowsGrid(),Grid.GetnColsGrid(),Grid.GetnMyRowsGrid(),Grid.GetnMyColsGrid());	
	}
	
	void DescInit(int *Des, int nM, int nN, CGrid& Grid)
	{
		Des[0] = 1;
		Des[1] = _nVal;
		Des[2] = nM;
		Des[3] = nN;
		Des[4] = Grid.GetnBlockSizeRows();
		Des[5] = Grid.GetnBlockSizeCols();
		Des[6] = 0;
		Des[7] = 0;
		Des[8] = Grid.GetnRowsDisGrid();		
	} 	
	
	int GetMyIdBLACS() const {return _nmyidBLACS;}
	int GetNumProcsBLACS() const {return _nnumprocsBLACS;}
	int GetnContinue() {return _nContinue;}
	int GetContxt(){return _nContxt;}
	int GetnWhat(){return _nWhat;}
	int GetnVal()  {return _nVal;}	
	
private:
	int _nmyidBLACS,_nnumprocsBLACS,_nContinue;
	int _nContxt,_nWhat, _nVal;
	char _order[3];	
};

class CExept
{
public:
	CExept(string pstring ):_string(pstring){}
	
	void traite()
	{
		std::cout << _string << std::endl;
	}
	
private:
	string _string;
};
