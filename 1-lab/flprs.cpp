#include <iostream>
#include <fstream>
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
double *distanceMatrix;	// Atstumu matrica.
int *X;					// Sprendinys

//=============================================================================

double getTime();
void loadDemandPoints();
void createDistanceMatrix();
void importMatrix();
void exportMatrix();
void randomSolution(int *X);
double HaversineDistance(double* a, double* b);
double evaluateSolution(int *X);

//=============================================================================

int main() {
	
	double t_start = getTime();          // Algoritmo vykdymo pradzios laikas

	loadDemandPoints();             // Nuskaitomi duomenys
	
	X = new int[numX];				// Sprendinys
	double u;						// Sprendinio tikslo funkcijos reiksme
	int *bestX = new int[numX];		// Geriausias rastas sprendinys
	double bestU = -1;				// Geriausio sprendinio tikslo funkcijos reiksme

	//----- Atstumu matricos sudarymas ----------------------------------------

	// Uncomment following two lines when running this code for the first time:
	createDistanceMatrix();
	exportMatrix();
	// importMatrix();
	double t_data_loaded = getTime();
	

	//----- Pagrindinis ciklas ------------------------------------------------
	
	for (int iters=0; iters<115000; iters++) {
		// Generuojam atsitiktini sprendini ir tikrinam ar jis nera geresnis uz geriausia zinoma
		randomSolution(X);
		u = evaluateSolution(X);
		if (u > bestU) {     // Jei geresnis, tai issaugojam kaip geriausia zinoma
			bestU = u;
			for (int i=0; i<numX; i++) bestX[i] = X[i];
			printf("Best result is: %f \n", bestU);
		}
		// if (iters % 100 == 0) printf("i = %d \n", iters);
	}
	//----- Rezultatu spausdinimas --------------------------------------------
	
	double t_calc_finished = getTime();     // Skaiciavimu pabaigos laikas

	cout << "Geriausias sprendinys: ";
	for (int i=0; i<numX; i++) cout << bestX[i] << " ";
	cout << "(" << bestU << ")" << endl << "Skaiciavimo trukme: " << t_calc_finished - t_data_loaded << endl;
	printf("Duomenu ikelimo trukme: %.5f \n", t_data_loaded - t_start);
	cout << "Sumine trukme: " <<  t_calc_finished - t_start << endl;

}

//=============================================================================

void loadDemandPoints() {
	
	//----- Load demand points ------------------------------------------------
	FILE *f;
	f = fopen("demandPoints.dat", "r");
	demandPoints = new double*[numDP];
	for (int i=0; i<numDP; i++) {
		demandPoints[i] = new double[3];
		fscanf(f, "%lf%lf%lf", &demandPoints[i][0], &demandPoints[i][1], &demandPoints[i][2]);
	}
	fclose(f);
}

//=============================================================================

void createDistanceMatrix() {
	distanceMatrix = new double[numDP*(numCL + 1)]; // new double[numDP*numDP];
	#pragma omp parallel for
	for (int i=0; i<numCL + 1; i++) 
		for (int j=0; j<numDP; j++) {
		distanceMatrix[i*numDP + j] = HaversineDistance(demandPoints[i], demandPoints[j]);	
	}
}

//=============================================================================

void importMatrix() {
	int size = numDP*(numCL + 1); // numDP*numDP;
	distanceMatrix = new double[size];
	ifstream rf("distance-matrix.dat", ios::out | ios::binary);
	if(!rf) {
      cout << "Failed to open distance matrix file." << endl;
      return;
   	}
	// Reading from the document in chunks of 10000.
	for (int i=0; i<numCL + 1; i++) {
		rf.read((char *) &distanceMatrix[i*numDP], sizeof(double[numDP]));
	}
	rf.close();
	
	if (!rf.good()) {
		cout << "Error occurred while reading distance matrix file." << endl;
	}
	printf("File distance-matrix was imported.\n");
}
//=============================================================================

void exportMatrix() {
	ofstream wf("distance-matrix.dat", ios::out | ios::binary);
	if(!wf) {
		cout << "Cannot open file!" << endl;
		return;
	}
	for (int i=0; i<numCL + 1; i++) {
		wf.write((char *) &distanceMatrix[i*numDP], sizeof(double[numDP]));
	}
	
	wf.close();

	if(!wf.good()) {
		cout << "Error occurred at writing time!" << endl;
	}
	printf("File distance-matrix was created.\n");
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
	// Iterating all demand points to calculate the sum of potential clients for the new object.
	omp_set_num_threads(4);
	#pragma omp parallel for reduction(+:U) private(d, bestPF, bestX)
	for (int i=0; i<numDP; i++) {
		// Smallest distance to Preexisting Facilities.
		bestPF = 1e5;
		for (int j=0; j<numPF; j++) {
			d = distanceMatrix[j*numDP + i]; // HaversineDistance(demandPoints[i], demandPoints[j]);
			if (d < bestPF) bestPF = d;
		}
		// Smallest distance to Candidate Locations.
		bestX = 1e5;
		for (int j=0; j<numX; j++) {
			d = distanceMatrix[X[j]*numDP + i]; // HaversineDistance(demandPoints[i], demandPoints[X[j]]);
			// cout << "Saved distance: " << distanceMatrix[X[j]*numDP + i] << ", calculated distance: " << HaversineDistance(demandPoints[i], demandPoints[X[j]]) << endl;
			if (d < bestX) bestX = d;
		}
		if (bestX < bestPF) {
			U += demandPoints[i][2];
		}
		else if (bestX == bestPF) U += 0.3*demandPoints[i][2];
	}
	// printf("Clients count: %f \n", U);
	return U;
}