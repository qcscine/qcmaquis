/* Fichier d'interface de la classe CRand Générateur de nombres aléatoires 
	donnés par COMBES */
using namespace std;

const int MAXI = 0x7FFF; 
class CRand
{	
	public:
	
	//Constructeur et destructeur
	
	CRand (int ninit=0);		//Le constructeur initialise un générateur
				// de nombres aléatoires	                			
	                        //Remplace randinit ()
	
	~CRand () {}	//Déclarée inline ne fait rien de spécial
	

	void InitRand (int nmyid)
	{
	  unsigned time_rand;
          long ltime;
	  time (&ltime);
	  time_rand = (unsigned int)(ltime&0x0000FFF);
          Randomize((nmyid+1)*time_rand);	
	}

	int Randint () 	//Retourne un entier entre 0 et MAXI
	{
		return nextrand (); 			
	}
	
	int Randn(int n) // Retourne un entier entre 0 et < n
	{
		return (nextrand()%n);
	}
	
	double Randdbl () //Retourne un nb aleatoire entre 0 et 1
	{
		int nVal = nextrand ();
		if (nVal == 0)
			return 1.0E-8;
		if (nVal == MAXI)
			return 0.9999999;
		return (double (nVal) / MAXI);
	}
	
	bool Randbool (double proba)
	{
		return ((nextrand () <= (int)(proba*MAXI)) ? true : false);
	}
	
	void Randomize (unsigned int n);
//	double Rand01 (int Idum);

	private:

	CRand (const CRand& Rand)
	{
	}
	
        CRand& operator = (const CRand& Rand)
        {
	 return *this;
        }

	int m_xi [17];
	
	int m_nind1,m_nind2;
	
	int nextrand () 		//Déclaré inline pour l'exécution rapide
	{
		m_nind1++;
		m_nind2++;
		m_nind1 %= 17;
		m_nind2 %= 17;
		return (m_xi [m_nind1] = ((m_xi [m_nind1] + m_xi [m_nind2]) & MAXI));
	}
	static int m_yi [17];
};	

const int Mbigcpp = 1000000000;
const int MSeedcpp = 161803398;
const int Mzcpp =0;
const double Faccpp = 1.0E-9;

class CRandB		//Cette classe CRand reprend les fonctionnalités de la classe CRand (Combes) avec le générateur de Bird
{
public:
	CRandB (int ninit=0);


	void init (unsigned int nInit);
	double Randdbl ()
	{
		double Rf;
		do
		{
			Rf = NextRand ()*Faccpp;
		}
		while ((Rf<= 1.0E-8) || (Rf >= 0.99999999));
		return Rf;
	}
	
	unsigned int RandInt (unsigned int nMax)
	{
		return NextRand () % nMax;
	}
	
	bool Randbool (double Proba)
	{
		double MaxProb = Mbigcpp*Proba;
		return((NextRand () <= MaxProb)?true:false);
	}
	
private:
       CRandB (const CRandB& RandB)
       {
       }

       CRandB& operator = (const CRandB& RandB )
       {
	 return *this;
       }


       int NextRand ()
	{
			int Mj;
			m_Inext++;
			m_Inextp++;
			m_Inext %= 55;
			m_Inextp %= 55;
			Mj = m_Ma [m_Inext] - m_Ma [m_Inextp];
			if (Mj < 0)
				Mj += Mbigcpp;
			m_Ma [m_Inext] = Mj;
			return Mj;
	}
	
	 int  m_Ma [55];
	 int  m_Inext;
	 int  m_Inextp;
};	

	
	
	
	
	
								
