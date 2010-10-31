#include <mpi.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <cstring>

#include "Crand.h"
#include "CMatrix.h"
#include "CPar.h"
#include "Const.h"

using namespace std;

CRandB* pRand;
/*
	Inside the main friendly user interface, C++ references
*/


int main (int argc, char *  argv[]) 
{	
	
/*
---------------------------- INIT MPI, BLACS
NOTE : all closed fonctions are done inside the two destructors 
*/

	int nM,nN,nSize;

	//size matrix, and of siwe of blocks 
	nM    = atoi(argv[1]);
	nN    = atoi(argv[2]);
	nSize = atoi(argv[3]);	

	// Init MPI & Blacs
	CParBLACS Par(argc,argv) ;
	// Init Grid
	CGrid Grid(Par.GetNumProcsMPI(),nSize);
	// Init coordinate of the grid, need MPI Parameters
	Par.InitGrid(Grid);
	// Init two cycl distrib
	Grid.InitBlockSize(nM,nN);
	
#ifdef DEBUG	
	Par.Print(Grid);
#endif
	
/*
---------------------------- START THE PDGESVD ROUTINE
*/
// Control parameters for blacs
	/*	
	WARNING Random generator is the same for every processor, if you need different, use next line
	srand((unsigned int) ((time ( NULL)+3*Par.GetNumProcsMPI())&0XFFF));	
*/
	srand((unsigned int) ((time ( NULL) &0XFFF)));
	time_t time_start, time_end;
	
	pRand = new CRandB(rand());
	
	Grid.InitBlockSize(nM,nN);

//Blacs control settings	
	int descA[9], descU[9], descV[9];

	Par.DescInit(descA,nM,nN,Grid);	
	Par.DescInit(descU,nM,nN,Grid);	
	Par.DescInit(descV,nM,nN,Grid);	

	CMatrix A(Grid.GetnRowsDisGrid(),Grid.GetnColsDisGrid());
	A.InitRandom();


	CMatrix U(Grid.GetnRowsDisGrid(),Grid.GetnColsDisGrid());
	CMatrix V(Grid.GetnRowsDisGrid(),Grid.GetnColsDisGrid());

	double *S = new double[max(nM,nN)];
	
	int IA=1,JA=1,IU=1,JU=1,IV=1,JV=1,info = 0;
	char JOBU[1] ={'V'};
	char JOBV[1] ={'V'};
	int lwork = -1;
	double *work = new double[1];
	double *rwork = new double[(1+4*max(nM,nN))];
	
	
	time_start = time(NULL);
	
	pdgesvd_(JOBU, JOBV, &nM, &nN,
		&A(0,0),&IA,&JA,descA,
		S,
		&U(0,0),&IU,&JU,descU,
		&V(0,0),&IV,&JV,descV,
		work, &lwork,
		rwork,
		&info);
	
	lwork = (int)(work[0]);
	
	work = new double[lwork];
	

	pdgesvd_(JOBU, JOBV, &nM, &nN,
		&A(0,0),&IA,&JA,descA,
		S,
		&U(0,0),&IU,&JU,descU,
		&V(0,0),&IV,&JV,descV,
		work, &lwork,
		rwork,
		&info);

	time_end = time(NULL);		

	if(Par.GetMyIdMPI() == 0)
	{
	ofstream out;
        out.open("time.txt",ios::app);
              out <<  time_end - time_start  <<endl;
        out.close();
	}

	
	delete[] S;
	delete[] work;
	delete[] rwork;


	return 0;
}
