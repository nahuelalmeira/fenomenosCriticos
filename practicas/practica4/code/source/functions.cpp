#include "../headers/functions.hpp"
#include "../headers/ran2.hpp"

#include <vector>

using std::vector;

double exponencial(vector<double> expVec, int dE)	//calcula exp(-dE/T)
{
	if(dE == 8) return expVec[1];
	if(dE == 4) return expVec[0];
	return 1;
}

void pasoMC(int L, vector<vector<int> > s, vector<double> expVec, int &E, int &m)
{
	int i, j;
	long p = 133;
	float eta;
	int dE;
	int h[L][L];
	int plus[L], minus[L];
	
	for(i = 0 ; i < L-1 ; i++)
		plus[i] = i+1;
	plus[L-1] = 0;
	for(i = 1 ; i < L ; i++)
		minus[i] = i-1;
	minus[0] = L-1;	
	
	for(i = 0 ; i < L ; i++)
	{
		for(j = 0 ; j < L ; j++)
		{
			h[i][j] = s[plus[i]][j] + s[minus[i]][j] + s[i][plus[j]] + s[i][minus[j]];
			dE = 2*s[i][j]*h[i][j];
			if(dE <= 0)
			{
				s[i][j] *= -1;
				m += 2*s[i][j];
				E += dE;
			}
			else
			{
				eta = ran2(&p);
				if(eta < exponencial(expVec, dE))
				{
					s[i][j] *= -1;
					m += 2*s[i][j];
					E += dE;
				}
			}
			
		}
	 }
}

