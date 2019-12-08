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
	int n_samples = atoi(argv[8]);


	printf("L   = %d\n", L);
	printf("L*L = %d\n", L*L);

	int i, j, k;    
	int samples;
	int **s;					//matriz de momentos magnéticos
	int plus[L], minus[L];
	double *expVec;				//almacena valores de exponenciales
	int E, M;
	double m;
	double T;
	double EMedia, mMedia, E2Media, m2Media, m4Media, cp, chi;
	FILE *gpipe;
	gpipe = popen("gnuplot --persist", "w");

	FILE*stream= fopen(argv[1], "w");
    assert(stream != NULL);
	
	expVec = AllocVector(2);
	s = AllocMatrixInt(L);
	E = 0;
	M = m = 0;
	
	for(i = 0 ; i < L-1 ; i++)
		plus[i] = i+1;
	plus[L-1] = 0;
	for(i = 1 ; i < L ; i++)
		minus[i] = i-1;
	minus[0] = L-1;	

	for(i=0; i<L; i++)
		for(j=0; j<L; j++)
			s[i][j] = 1;
	
	//Energía inicial
	for(i = 0 ; i < L ; i++)
		for(j = 0 ; j < L ; j++)
			E = E - s[i][j]*s[plus[i]][j] - s[i][j]*s[i][plus[j]];
	printf("E = %d\n", E);
	
	//Magnetización inicial
	for(i = 0 ; i < L ; i++)
		for(j = 0 ; j < L ; j++)
			M += s[i][j];
	printf("M = %d\n", M);
			
	fprintf(gpipe, "plot '-' w l lw 2\n");
	
	printf("T           e        m          m2          m4\n");

	//Barrido en temperatura
	for(T = Ti ; T <= Tf ; T += deltaT)
	{
		
		expVec[0] = exp(-4./T);
		expVec[1] = exp(-8./T);
		for(k = 0 ; k <= N ; k++) 
			pasoMC(L, s, expVec, &E, &M);
	
		EMedia  = 0;
		E2Media = 0;
		mMedia  = 0;
		m2Media = 0;
		m4Media = 0;
		samples = 0;
		
		for(k=0; k<Ns; k++)
		{
			pasoMC(L, s, expVec, &E, &M);
			if (k%(Ns/n_samples) == 0) {
				m = (1.*abs(M))/(L*L);
				EMedia += E;
				E2Media += E*E;
				mMedia += m;
				m2Media += m*m;
				m4Media += pow(m, 4);
				//printf("    %.3lf    %.8lf     %.8lf\n", mMedia, m2Media, m4Media);
				samples++;
			}
		}
		assert(samples==n_samples);
		
		//Divide por número de pasos MC
		EMedia  /= n_samples;
		E2Media /= n_samples;
		mMedia  /= n_samples;
		m2Media /= n_samples;
		m4Media /= n_samples;
		
		//Calcula calor específico y susceptibilidad magnética
		cp = (E2Media - EMedia*EMedia)/(L*L*T*T);
		chi = L*L*(m2Media - mMedia*mMedia)/T;	

		fprintf(gpipe, "%g %g\n", T, mMedia);
		//fprintf(stream, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",T,EMedia, EMedia/(L*L), mMedia, mMedia/(L*L), cp, chi);	
		fprintf(stream, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
				T, EMedia/(L*L), mMedia, cp, chi, m2Media, m4Media);	

		printf("%.3lf    %.3lf     %.3lf    %.8lf     %.8lf\n", T, EMedia/(L*L), mMedia, m2Media, m4Media);
		
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
