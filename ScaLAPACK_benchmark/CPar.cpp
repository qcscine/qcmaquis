#include "CPar.h"

void InitGrid::InitGrid(CGrid& Grid)
{
	Cblacs_gridinit(GetnContxt(),'R',Grid.GetnRows(),Grid.GetnCols());
};
