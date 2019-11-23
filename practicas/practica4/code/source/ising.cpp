#include "../headers/functions.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using std::cout;
using std::endl;
using std::vector;
using std::fill;
using std::ofstream;

//define N 100000		//Iteraciones MC hasta llegar al equilibrio
//define Ns 100000		//Iteraciones MC para calcular valores medios

int main(int argc, char *argv[])
{

	// Command line arguments
	//char output_file = argv[1];
	int L = atoi(argv[2]);
	int N = atoi(argv[3]);
	int Ns = atoi(argv[4]);
	float Ti = atof(argv[5]);
	float Tf = atof(argv[6]);
	float deltaT = atof(argv[7]);

	cout << "L = " << L << endl;

	int i, j, k;    
	vector<vector<int> > s(L, vector<int>(L, 1));				//matriz de momentos magnéticos
    int plus[L], minus[L];
	vector<double> expVec;				//almacena valores de exponenciales
	int E, m;
	double T;
	double EMedia, mMedia, E2Media, m2Media, m4Media, cp, chi;
	

	ofstream output;
    output.open(argv[1]);
    //assert(stream != NULL);
	
	E = 0;
	m = 0;
	
	for(i = 0 ; i < L-1 ; i++)
		plus[i] = i+1;
	plus[L-1] = 0;
	for(i = 1 ; i < L ; i++)
		minus[i] = i-1;
	minus[0] = L-1;	
	
	//Energía inicial
	for(i = 0 ; i < L ; i++)
		for(j = 0 ; j < L ; j++)
			E = E - s[i][j]*s[plus[i]][j] - s[i][j]*s[i][plus[j]];
	printf("E = %d\n", E);
	
	//Magnetización inicial
	for(i = 0 ; i < L ; i++)
		for(j = 0 ; j < L ; j++)
			m += s[i][j];
	
	//Barrido en temperatura
	for(T = Ti ; T <= Tf ; T += deltaT)
	{
		printf("Temperatura = %.2lf\n", T);
		expVec.push_back(exp(-4./T));
		expVec.push_back(exp(-8./T));
		for(k = 0 ; k <= N ; k++) 
			pasoMC(L, s, expVec, E, m);
	
		EMedia = 0;
		E2Media = 0;
		mMedia = 0;
		m2Media = 0;
		m4Media = 0;
		
		for(k = 0 ; k < Ns ; k++)
		{
			pasoMC(L, s, expVec, E, m);
			EMedia += E;
			E2Media += E*E;
			mMedia += abs(m);
			m2Media += m*m;
			m4Media += m*m*m*m;
		}
		
		//Divide por número de pasos MC
		EMedia /= Ns;
		E2Media /= Ns;
		mMedia = mMedia / Ns;
		m2Media /= Ns;
		m4Media /= Ns;
		
		//Calcula calor específico y susceptibilidad magnética
		cp = (E2Media - EMedia*EMedia)/(L*L*T*T);
		chi = (m2Media - mMedia*mMedia)/(L*L*T);	

		//fprintf(stream, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",T,EMedia, EMedia/(L*L), mMedia, mMedia/(L*L), cp, chi);	
		output << T  << "\t" << EMedia/(L*L) << "\t" << mMedia/(L*L)     << "\t" 
               << cp << "\t" << chi          << "\t" << m2Media/pow(L,4) << "\t" << m4Media/pow(L,8) 
               << endl;	
		
	}	
    output.close();
	return 0;
}