#include <iostream>
#include <vector>
#include <ctime>
#include <fstream>
#include <string>

#include "Const.h"
#include "CSMatrix.cpp"

using namespace std;


int main (int argc, char * const argv[]) 
{  
	//Datas
	int	NUM_ROWS(static_cast<size_t> (atoi(argv[1])));
	int	NUM_COLS(static_cast<size_t> (atoi(argv[2])));	
	string file;
	string extention(".txt");
	
	time_t time_start, time_end;
	ofstream out;
	
/*------------------------------------------------ DGEMM Benchmark--------------------------------------------------------------*/	

#ifdef DGEMM
	{ //<- tip for the memory
	//outputfile
	file = "time_DGEMM_";	
	file = file + argv[4] + extention;
	
	//Matrix NUM_ROWS >= NUM_COLS
	CSMatrix<double> A(NUM_ROWS,NUM_COLS);
	CSMatrix<double> B(NUM_COLS,NUM_ROWS);
	CSMatrix<double> C(NUM_ROWS,NUM_ROWS);

	double one = 1;
	double zero = 0;
	
	time_start = time(NULL);

	//dgemm_("T",  "T", &A.NUM_ROWS(), &NUM_ROWS, &NUM_ROWS, &one, &A(0,0), &NUM_ROWS, &B(0,0), &NUM_COLS, &zero, &C(0,0), &NUM_ROWS);
		
	int ANumRow = A.GetnNumRow();
	int ANumCol = A.GetnNumCol();
	int BNumRow = B.GetnNumRow();		
	int CNumRow = C.GetnNumRow();
		
		
	dgemm_("T",  "T", &ANumRow, &BNumRow, &ANumCol, &one, &A(0,0), &ANumRow, &B(0,0), &BNumRow, &zero, &C(0,0), &CNumRow);

	time_end = time(NULL);
	
	
	out.open(file.c_str(),ios::app);
		out <<  time_end - time_start << " " << argv[3]  << " " << argv[4] << endl;
	out.close();
	}
#endif	

/*------------------------------------------------ DGESVD Benchmark--------------------------------------------------------------*/	
	
#ifdef DGESVD	
	{ //<- tip for the memory
	//outputfile
	file = "time_DGESVD_";	
	file = file + argv[4] + extention;	
		
	//Matrix
	CSMatrix<double> A(NUM_ROWS,NUM_COLS);
	CSMatrix<double> U(NUM_ROWS,NUM_ROWS);
	CSMatrix<double> VT(NUM_COLS,NUM_COLS);	
	CSVector<double> S(NUM_COLS);
	
	int lwork = -1;
	int info = 0;
	double wkopt;
	
	int ANumRow = A.GetnNumRow();
	int ANumCol = A.GetnNumCol();
	int UNumRow = U.GetnNumRow();		
	int VTNumRow = VT.GetnNumRow();					
		
	time_start = time(NULL);	
		
	dgesvd_("A", "A",&ANumRow,&ANumCol,&A(0,0),&ANumRow,&S(0),&U(0,0),&UNumRow,&VT(0,0),&VTNumRow,&wkopt,&lwork,&info); 
		
	lwork = static_cast<long int> (wkopt);
	CSVector<double> work(lwork);
	
	dgesvd_("A", "A",&ANumRow,&ANumCol,&A(0,0),&ANumRow,&S(0),&U(0,0),&UNumRow,&VT(0,0),&VTNumRow,&work(0),&lwork,&info); 
			
	time_end = time(NULL);			
		
	out.open(file.c_str(),ios::app);
		out <<  time_end - time_start << " " << argv[3]  << " " << argv[4] << endl;
	out.close();					
	}	
#endif	
	
/*------------------------------------------------ DGESDD Benchmark--------------------------------------------------------------*/	
	
#ifdef DGESDD	
	{ //<- tip for the memory
		//outputfile
		file = "time_DGESDD_";	
		file = file + argv[4] + extention;	
		
		//Matrix
		CSMatrix<double> A(NUM_ROWS,NUM_COLS);
		CSMatrix<double> U(NUM_ROWS,NUM_ROWS);
		CSMatrix<double> VT(NUM_COLS,NUM_COLS);	
		CSVector<double> S(NUM_COLS);
		CSVector<int> IWORK(8*min(NUM_ROWS,NUM_COLS));

		int lwork = -1;
		int info = 0;
		double wkopt;
		
		int ANumRow = A.GetnNumRow();
		int ANumCol = A.GetnNumCol();
		int UNumRow = U.GetnNumRow();		
		int VTNumRow = VT.GetnNumRow();					
		
		time_start = time(NULL);	

		dgesdd_("A",&ANumRow,&ANumCol,&A(0,0),&ANumRow,&S(0),&U(0,0),&UNumRow,&VT(0,0),&VTNumRow,&wkopt,&lwork,&IWORK(0),&info);
		
		lwork = static_cast<long int> (wkopt);
		CSVector<double> work(lwork);
		
		dgesdd_("A",&ANumRow,&ANumCol,&A(0,0),&ANumRow,&S(0),&U(0,0),&UNumRow,&VT(0,0),&VTNumRow,&work(0),&lwork,&IWORK(0),&info);

		time_end = time(NULL);			
		
		out.open(file.c_str(),ios::app);
			out <<  time_end - time_start << " " << argv[3]  << " " << argv[4] << endl;
		out.close();					
	}	
#endif	
	
	
    return 0;
}
