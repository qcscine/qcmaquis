#include <iostream>
#include <vector>
#include <ctime>


#include "CSMatrix.h"

using namespace std;



int main (int argc, char * const argv[]) 
{  
	int n=1;
	size_t NUM(40000);
	
	CSMatrix A;

	{

		time_t tstart, tend;
	
		
		CSMatrix B(NUM,NUM);
		B.Init(0);		
		A=B;
		tstart = time(NULL);
		A(1,0) = 10;
		tend = time(NULL);
		
		cout << tend - tstart << endl;
	}
	
	//A.Print();

/*
	cout << "Vectors : " << endl;	
	
	{
		CSVector A(2);
		
		{
			CSVector B(2);
			B.Init(1);		
			B.Print();
			A=B;
			A.Print();
			A(1) = 10;
		//	A.Print();		
		}
		
		A.Print();
	}
*/	
	
	
    return 0;
}
