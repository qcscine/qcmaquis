#include <iostream>
#include <cmath>
#include <ctime>
#include <cstring>

#include "CMatrix.h"
#include "Crand.h"
#include "Const.h"


void CMatrix::InitRandom()
{
	for(int i=0;i<_nColumns;i++)
	{
		for(int j=0;j<_nRows;j++)
		{
			_ppM[i][j] = pRand->Randdbl();
		}
	}
}

