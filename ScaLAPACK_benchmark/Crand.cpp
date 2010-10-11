//Implémentation de la classe CRand fichier CRAND.CPP

#include <iostream>
#include <cmath>

#include "Crand.h"

using namespace std;


int CRand::m_yi [17] = {22934,3922,10976,32327,13781,9407,20980,32338,6184,
							19583,1676,6438,31235,29580,13181,20665,20441};

	
CRand::CRand(int ninit)
         :m_nind1 (16),m_nind2 (11)
{

	for (register int i = 0; i<17; i++)
		m_xi [i] = m_yi [i];
	
        Randomize (ninit);		
}

void CRand::Randomize (unsigned int n)
{
	unsigned i;
	unsigned n1 = n&0xFFFF;
	for (i = 0; i<n1 ; i++)
		nextrand ();
}			


CRandB::CRandB (int ninit)
{	
	for (int i = 0;i < 55;i++)
		m_Ma [i] = 0;
	m_Inext = 0;
	m_Inextp = 31;
	init(ninit);
}

void CRandB::init (unsigned int nInit)
{
	int Mj,Mk,Ii;
	Mj = MSeedcpp - nInit;
	Mj %= Mbigcpp;
	m_Ma [54] = Mj;
	Mk = 1;
	for (int i = 1;i<55;i++)
	{
		Ii = 21*i%55;
		m_Ma [Ii-1] = Mk;
		Mk = Mj - Mk;
		if (Mk < Mzcpp)
			Mk += Mbigcpp;
		Mj = m_Ma [Ii-1];
	}
	for (int k = 1;k<=4;k++)
	{
		for (int i = 0; i<55;i++)
		{
			m_Ma [i] -= m_Ma [(i+30)%55];
			if (m_Ma[i] < Mzcpp)
				m_Ma [i] += Mbigcpp;
		}
	}
	m_Inext = 0;
	m_Inextp = 31;
	NextRand ();
}
