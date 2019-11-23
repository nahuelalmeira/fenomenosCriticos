#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <complex.h>
#include <stdio.h>

//define N 100000		//Iteraciones MC hasta llegar al equilibrio
//define Ns 100000		//Iteraciones MC para calcular valores medios

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran2(long *);
double exponencial(double[], int);	//calcula exp(-dE/T)
double *AllocVector(int);			//alloca vector double
int **AllocMatrixInt(int);			//alloca matriz entero
void FreeMatrixInt(int **, int);	//libera matriz entero
void pasoMC(int L, int **s, double *expVec, int *E, int *m);

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

	printf("L = %d\n", L);

	int i, j, k;    
	int **s;					//matriz de momentos magnéticos
	int plus[L], minus[L];
	double *expVec;				//almacena valores de exponenciales
	int E, m;
	double T;
	double EMedia, mMedia, E2Media, m2Media, m4Media, cp, chi;
	FILE *gpipe;
	gpipe = popen("gnuplot --persist", "w");

	FILE*stream= fopen(argv[1], "w");
    assert(stream != NULL);
	
	expVec = AllocVector(2);
	s = AllocMatrixInt(L);
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
			
	fprintf(gpipe, "plot '-' w l lw 2\n");
	
	//Barrido en temperatura
	for(T = Ti ; T <= Tf ; T += deltaT)
	{
		printf("Temperatura = %.2lf\n", T);
		expVec[0] = exp(-4./T);
		expVec[1] = exp(-8./T);
		for(k = 0 ; k <= N ; k++) 
			pasoMC(L, s, expVec, &E, &m);
	
		EMedia = 0;
		E2Media = 0;
		mMedia = 0;
		m2Media = 0;
		m4Media = 0;
		
		for(k = 0 ; k < Ns ; k++)
		{
			pasoMC(L, s, expVec, &E, &m);
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

		fprintf(gpipe, "%g %g\n",T,mMedia/(L*L));
		//fprintf(stream, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",T,EMedia, EMedia/(L*L), mMedia, mMedia/(L*L), cp, chi);	
		fprintf(stream, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
				T, EMedia/(L*L), mMedia/(L*L), cp, chi, m2Media/pow(L,4), m4Media/pow(L,8));	
		
	}
	fprintf(gpipe, "e\n");	
	FreeMatrixInt(s, L);
	free(expVec);	
	return 0;
}

double *AllocVector(int n)
{
    double *v = (double *) malloc(n * sizeof(double));
    assert(v != NULL);
    return v;
}

int **AllocMatrixInt(int n)
{
    int i;

    int **m = (int **) malloc(n * sizeof(int*));
    assert(m != NULL);
    for( i = 0; i < n; i++ )
        m[i] = (int *) malloc(n * sizeof(int));
    return m;
}

void FreeMatrixInt(int **pMat, int n)
{
    int i;

    for( i = 0; i < n; i++ )
        free(pMat[i]);
    free(pMat);
    return;
}

double exponencial(double expVec[2], int dE)	//calcula exp(-dE/T)
{
	if(dE == 8) return expVec[1];
	if(dE == 4) return expVec[0];
	return 1;
}

void pasoMC(int L, int **s, double *expVec, int *E, int *m)
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
				*m += 2*s[i][j];
				*E += dE;
			}
			else
			{
				eta = ran2(&p);
				if(eta < exponencial(expVec, dE))
				{
					s[i][j] *= -1;
					*m += 2*s[i][j];
					*E += dE;
				}
			}
			
		}
	 }
}

float ran2(long *idum){
int j;
long k;
static long idum2=123456789;
static long iy=0;
static long iv[NTAB];
float temp;
if (*idum <= 0) {
if (-(*idum) < 1) *idum=1;
else *idum = -(*idum);
idum2=(*idum);
for (j=NTAB+7;j>=0;j--) { 
k=(*idum)/IQ1;
*idum=IA1*(*idum-k*IQ1)-k*IR1;
if (*idum < 0) *idum += IM1;
if (j < NTAB) iv[j] = *idum;
}
iy=iv[0];
}
k=(*idum)/IQ1;
*idum=IA1*(*idum-k*IQ1)-k*IR1;
if (*idum < 0) *idum += IM1;
k=idum2/IQ2;
idum2=IA2*(idum2-k*IQ2)-k*IR2;
if (idum2 < 0) idum2 += IM2;
j=iy/NDIV;
iy=iv[j]-idum2;
iv[j] = *idum;
if (iy < 1) iy += IMM1;
if ((temp=AM*iy) > RNMX) return RNMX;
else return temp;
}
