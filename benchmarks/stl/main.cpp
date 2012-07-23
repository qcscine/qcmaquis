#include <iostream>
#include <vector>
#include <time.h>

#include "cow_vector.h"
#include "aligned_allocator.h"

#include "CSMatrix.h"

using namespace std;



int main (int argc, char * const argv[]) 
{  
	int n=1;
	size_t NUM(33333);
	

/*
	{
		
		CSMatrix A;

		time_t tstart, tend;
	
		CSMatrix B(NUM,NUM);
		B.Init(0);		
		A=B;
		tstart = time(NULL);
		A(0,0) = 10;
		tend = time(NULL);
		
		cout << " Tim :" << tend - tstart << endl;
	}
	
	
	{
		
		time_t tstart, tend;
		
		//copy_on_write_vector<double, aligned_allocator<double, 8> > B(NUM*NUM);
		//copy_on_write_vector<double, aligned_allocator<double, 8, true> > A;
		
		copy_on_write_vector<double> B(NUM*NUM), A;
		
		A=B;
		tstart = time(NULL);
		A[0] = 10;
		tend = time(NULL);
		
		cout << " Bella :" << tend - tstart << endl;
		
		
	}
*/
	
	{
		clock_t tstart, tend;
		tstart = clock();
				//copy_on_write_vector<double, aligned_allocator<double, 8> > B(NUM*NUM);
				vector<double> B(NUM*NUM);
		tend = clock();
		
		cout << " Bella 0 :" << (tend - tstart) << "micro s" <<  endl;

		tstart = clock();
			vector<double, aligned_allocator<double, 8, true> > A(NUM*NUM);
		tend = clock();
		
		cout << " Bella 1 :" << (tend - tstart) << "micro s" << endl;
	}
	
    return 0;
}
