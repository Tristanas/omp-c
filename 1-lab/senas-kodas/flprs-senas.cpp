#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>

using namespace std;

int numDP = 10000;      // Vietoviu skaicius (demand points, max 10000)
int numPF = 5;          // Esanciu objektu skaicius (preexisting facilities)
int numCL = 50;         // Kandidatu naujiems objektams skaicius (candidate locations)
int numX  = 3;          // Nauju objektu skaicius

double **demandPoints;  // Geografiniai duomenys
int *X;					// Sprendinys
double **pfDistances;

double genFileTime = 0;
double openFileTime = 0;

//=============================================================================

void GenDistanceFile();
void LoadDistanceFile();
double getTime();
void loadDemandPoints(); // Uzkrauna dokumenta is aukstesnio folderio.
void randomSolution(int *X);
double HaversineDistance(double* a, double* b);
double evaluateSolution(int *X);

//=============================================================================

int main() {
	
	double ts = getTime();          // Algoritmo vykdymo pradzios laikas

	loadDemandPoints();             // Nuskaitomi duomenys

	// GenDistanceFile();            // Uzkomentuoti po pirmo nuskaitymo
	LoadDistanceFile();

	/*pfDistances = new double*[numCL];
	for (int i = 0; i < numCL; i++)
	{
		pfDistances[i] = new double[numDP];
	}
	
	//#pragma omp parallel for schedule(guided) collapse(2)
	for(int i = 0; i < numCL;i++)
	{
		
		for(int j = 0; j < numDP; j++)
			pfDistances[i][j] = HaversineDistance(demandPoints[j], demandPoints[i]);
	}*/

	double tska = getTime();

	int nProcessors = 4;
	omp_set_num_threads(nProcessors);

	
	double u;						// Sprendinio tikslo funkcijos reiksme
	int *bestX = new int[numX];		// Geriausias rastas sprendinys
	double bestU = -1;				// Geriausio sprendinio tikslo funkcijos reiksme

	
	//----- Pagrindinis ciklas ------------------------------------------------
	#pragma omp parallel for schedule(guided) private(u, X)
	for (int iters=0; iters<2848; iters++) {
		// Generuojam atsitiktini sprendini ir tikrinam ar jis nera geresnis uz geriausia zinoma
		X = new int[numX];
		randomSolution(X);
		u = evaluateSolution(X);
		if (u > bestU) {     // Jei geresnis, tai issaugojam kaip geriausia zinoma
			bestU = u;
			for (int i=0; i<numX; i++) bestX[i] = X[i];
		}
	}
	//----- Rezultatu spausdinimas --------------------------------------------
	
	double tf = getTime();     // Skaiciavimu pabaigos laikas

	cout << "Geriausias sprendinys: ";
	for (int i=0; i<numX; i++) cout << bestX[i] << " ";
	cout << "(" << bestU << ")" << endl 
	<< "Programos trukme: " << tf-ts << endl 
	<< "Skaiciavimo trukme: " << tf - tska << endl 
	<< "Failo generavimo trukme: " << genFileTime << endl 
	<< "Failo krovimo trukme: " << openFileTime << endl;
}

//=============================================================================

void GenDistanceFile()
{
	double ts = getTime();
	FILE *f;
	f = fopen("distances.dat", "w+");
	for(int i = 0; i < numCL;i++)
	{
		for(int j = 0; j < numDP; j++)
		{
			fprintf(f, "%.3f ", HaversineDistance(demandPoints[j], demandPoints[i]));
		}
	}
	fclose(f);
	genFileTime = getTime() - ts;
}

void LoadDistanceFile()
{
	int arraySize = numCL;
	double ts = getTime();
	FILE *f;
	f = fopen("distances.dat", "r");
	pfDistances = new double*[arraySize];
	for(int i = 0; i < arraySize;i++)
	{
		pfDistances[i] = new double[numDP];
		for(int j = 0; j < numDP; j++)
		{
			fscanf(f, "%lf", &pfDistances[i][j]);
		}
	}
	fclose(f);
	openFileTime = getTime() - ts;
}


//=============================================================================

void loadDemandPoints() {
	
	//----- Load demand points ------------------------------------------------
	FILE *f;
	f = fopen("../demandPoints.dat", "r");
	demandPoints = new double*[numDP];
	for (int i=0; i<numDP; i++) {
		demandPoints[i] = new double[3];
		fscanf(f, "%lf%lf%lf", &demandPoints[i][0], &demandPoints[i][1], &demandPoints[i][2]);
	}
	fclose(f);
}

//=============================================================================

double HaversineDistance(double* a, double* b) {
   double dlon = fabs(a[0] - b[0]);
   double dlat = fabs(a[1] - b[1]);
   double aa = pow((sin((double)dlon/(double)2*0.01745)),2) + cos(a[0]*0.01745) * cos(b[0]*0.01745) * pow((sin((double)dlat/(double)2*0.01745)),2);
   double c = 2 * atan2(sqrt(aa), sqrt(1-aa));
   double d = 6371 * c; 
   return d;
}

//=============================================================================

double getTime() {
   struct timeval laikas;
   gettimeofday(&laikas, NULL);
   double rez = (double)laikas.tv_sec+(double)laikas.tv_usec/1000000;
   return rez;
}

//=============================================================================

void randomSolution(int *X) {
	int unique;
	for (int i=0; i<numX; i++) {
		do {
			unique = 1;
			X[i] = (int)((double)rand()/RAND_MAX * numCL);
			for (int j=0; j<i; j++)
				if (X[j] == X[i]) {
					unique = 0;
					break;
				}		
		} while (unique == 0);
	}
}

//=============================================================================

double evaluateSolution(int *X) {
	double U = 0;
	int bestPF;
	int bestX;
	double d;
	for (int i=0; i<numDP; i++) {
		bestPF = 1e5;
		for (int j=0; j<numPF; j++) {
			d = pfDistances[j][i];//HaversineDistance(demandPoints[i], demandPoints[j]);
			if (d < bestPF) bestPF = d;
		}
		bestX = 1e5;
		for (int j=0; j<numX; j++) {
			d = pfDistances[X[j]][i];//HaversineDistance(demandPoints[i], demandPoints[X[j]]);
			if (d < bestX) bestX = d;
		}
		if (bestX < bestPF) U += demandPoints[i][2];
		else if (bestX == bestPF) U += 0.3*demandPoints[i][2];
	}
	return U;
}