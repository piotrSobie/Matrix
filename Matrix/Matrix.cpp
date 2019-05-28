#include "pch.h"
#include <iostream>
#include <cmath>
#include <ctime>
#include <fstream>

#define INDEKS 111111

using namespace std;

struct MatrixStruct {
	int m, n;
	double **matrix;
	MatrixStruct(int a, int b) {
		m = a;//m wiersze
		n = b;//n kolumny
		matrix = new double*[m];
		for (int i = 0; i < m; i++) {
			matrix[i] = new double[n];
		}
	}
};

void deleteMatrix(MatrixStruct *X) {
	for (int i = 0; i < X->m; i++) {
		delete[]X->matrix[i];
	}
	delete[]X->matrix;
}

void printMatrix(MatrixStruct *X) {
	if (X == NULL) {
		cout << "Matrix is null pointer" << endl;
		return;
	}
	for (int x = 0; x < X->m; x++) {
		for (int y = 0; y < X->n; y++) {
			cout << X->matrix[x][y] << ' ';
		}
		cout << endl;
	}
}

void fillMatrix(MatrixStruct *X, double a1, double a2, double a3, double rest) {
	for (int x = 0; x < X->m; x++) {
		for (int y = 0; y < X->n; y++) {
			if (x == y) {
				X->matrix[x][y] = a1;
			}
			else if (x == y - 1 || y == x - 1) {
				X->matrix[x][y] = a2;
			}
			else if (x == y - 2 || y == x - 2) {
				X->matrix[x][y] = a3;
			}
			else {
				X->matrix[x][y] = rest;
			}
		}
	}
}

void fillMatrix(MatrixStruct *X, double a1) {
	for (int x = 0; x < X->m; x++) {
		for (int y = 0; y < X->n; y++) {
			X->matrix[x][y] = a1;
		}
	}
}

MatrixStruct* addMatrix(MatrixStruct *A, MatrixStruct *B) {
	if (A->m != B->m || A->n != B->n) {
		return NULL;
	}
	MatrixStruct *Result = new MatrixStruct(A->m, A->n);
	for (int x = 0; x < Result->m; x++) {
		for (int y = 0; y < Result->n; y++) {
			Result->matrix[x][y] = A->matrix[x][y] + B->matrix[x][y];
		}
	}
	return Result;
}

MatrixStruct* minusMatrix(MatrixStruct *A, MatrixStruct *B) {
	if (A->m != B->m || A->n != B->n) {
		return NULL;
	}
	MatrixStruct *Result = new MatrixStruct(A->m, A->n);
	for (int x = 0; x < Result->m; x++) {
		for (int y = 0; y < Result->n; y++) {
			Result->matrix[x][y] = A->matrix[x][y] - B->matrix[x][y];
		}
	}
	return Result;
}

MatrixStruct* multiMatrix(MatrixStruct *A, MatrixStruct *B) {
	if (A->n != B->m) {
		return NULL;
	}
	MatrixStruct *Result = new MatrixStruct(A->m, B->n);
	for (int x = 0; x < Result->m; x++) {
		for (int y = 0; y < Result->n; y++) {
			Result->matrix[x][y] = 0;
			for (int z = 0; z < Result->m; z++) {
				Result->matrix[x][y] += A->matrix[x][z] * B->matrix[z][y];
			}
		}
	}
	return Result;
}

MatrixStruct* numberMultiMatrix(MatrixStruct *A, double a) {
	MatrixStruct *Result = new MatrixStruct(A->m, A->n);
	for (int x = 0; x < A->m; x++) {
		for (int y = 0; y < A->n; y++) {
			if (A->matrix[x][y] != 0) {
				Result->matrix[x][y] = A->matrix[x][y] * a;
			}
			else {
				Result->matrix[x][y] = A->matrix[x][y];
			}
		}
	}
	return Result;
}

double calculateNorm(MatrixStruct *A, MatrixStruct *x, MatrixStruct *b) {
	MatrixStruct *AxX = multiMatrix(A, x);
	MatrixStruct *res = minusMatrix(AxX, b);
	deleteMatrix(AxX);
	delete AxX;
	double norm = 0;
	for (int x = 0; x < res->m; x++) {
		for (int y = 0; y < res->n; y++) {
			norm += pow(res->matrix[x][y], 2);
		}
	}
	norm = sqrt(norm);

	deleteMatrix(res);
	delete res;

	return norm;
}

void fillMatrixULD(MatrixStruct *A, MatrixStruct *U, MatrixStruct *L, MatrixStruct *D) {
	for (int x = 0; x < A->m; x++) {
		for (int y = 0; y < A->n; y++) {
			if (x == y) {
				D->matrix[x][y] = A->matrix[x][y];
			}
			else {
				D->matrix[x][y] = 0;
			}

			if (x > y) {
				L->matrix[x][y] = A->matrix[x][y];
			}
			else {
				L->matrix[x][y] = 0;
			}

			if (y > x) {
				U->matrix[x][y] = A->matrix[x][y];
			}
			else {
				U->matrix[x][y] = 0;
			}
		}
	}
}

//Jacobi
//Ax=b 
//A=L+U+D
//(L+U+D)x=b
//Dx=-(L+U)x+b
//Dx=y	y=-(L+U)x+b
void Jacob(MatrixStruct *A, MatrixStruct *x, MatrixStruct *b, fstream& file) {
	MatrixStruct *D = new MatrixStruct(A->m, A->n); //macierz diagonalna
	MatrixStruct *L = new MatrixStruct(A->m, A->n); //macierz trojkatna dolna
	MatrixStruct *U = new MatrixStruct(A->m, A->n); //macierz trojkatna gorna

	fillMatrixULD(A, U, L, D);

	cout << "Jacobi:" << endl;
	file << "Jacobi:" << endl;
	int iJ = 0;
	clock_t start = clock();
	double norm = calculateNorm(A, x, b);
	double lastNorm;
	while (norm > pow(10, -9)) {
		lastNorm = norm;
		MatrixStruct *LpU = addMatrix(L, U);
		MatrixStruct *mLpU = numberMultiMatrix(LpU, -1);
		deleteMatrix(LpU);
		delete LpU;
		MatrixStruct *mLpUxX = multiMatrix(mLpU, x);
		deleteMatrix(mLpU);
		delete mLpU;
		MatrixStruct *y = addMatrix(mLpUxX, b);
		deleteMatrix(mLpUxX);
		delete mLpUxX;
		for (int xx = 0; xx < D->m; xx++) {
			for (int yy = 0; yy < D->n; yy++) {
				if (xx == yy) {
					x->matrix[xx][0] = y->matrix[xx][0] / D->matrix[xx][yy];
				}
			}
		}
		iJ++;
		deleteMatrix(y);
		delete y;
		norm = calculateNorm(A, x, b);
		if (norm > lastNorm) {
			cout << "Metoda Jacobiego nie zbiega sie" << endl;
			file << "Metoda Jacobiego nie zbiega sie" << endl;
			return;
		}
	}

	deleteMatrix(D);
	delete D;
	deleteMatrix(L);
	delete L;
	deleteMatrix(U);
	delete U;

	double finalNorm = calculateNorm(A, x, b);
	long time = (clock() - start);

	cout << "Ilosc iteracji: " << iJ << endl;
	cout << "Norma: " << finalNorm << endl;
	cout << "Czas: " << time << " ms." << endl;

	file << "Ilosc iteracji: " << iJ << endl;
	file << "Norma: " << finalNorm << endl;
	file << "Czas: " << time << " ms." << endl;
}

//Gauss-Seidel
//Ax=b
//A=L+U+D
//(L+U+D)x=b
//(D+L)x=-Ux+b
//y=-Ux+b	Z=D+L
//Zx=y
void GaussSeidel(MatrixStruct *A, MatrixStruct *x, MatrixStruct *b, fstream& file) {
	MatrixStruct *D = new MatrixStruct(A->m, A->n); //macierz diagonalna
	MatrixStruct *L = new MatrixStruct(A->m, A->n); //macierz trojkatna dolna
	MatrixStruct *U = new MatrixStruct(A->m, A->n); //macierz trojkatna gorna

	fillMatrixULD(A, U, L, D);

	cout << "Gauss-Seidel:" << endl;
	file << "Gauss-Seidel:" << endl;
	int iGS = 0;
	clock_t start = clock();
	double norm = calculateNorm(A, x, b);
	double lastNorm;
	double help;
	while (norm > pow(10, -9)) {
		lastNorm = norm;
		//cout << iGS << '\t' << calculateNorm(A, x, b) << endl;
		MatrixStruct *mU = numberMultiMatrix(U, -1);
		MatrixStruct *mUxX = multiMatrix(mU, x);
		deleteMatrix(mU);
		delete mU;
		MatrixStruct *y = addMatrix(mUxX, b);
		deleteMatrix(mUxX);
		delete mUxX;
		MatrixStruct *Z = addMatrix(D, L);
		for (int xx = 0; xx < Z->m; xx++) {
			help = y->matrix[xx][0];
			for (int zz = 0; zz < xx; zz++) {
				help -= (Z->matrix[xx][zz]) * (x->matrix[zz][0]);
			}
			x->matrix[xx][0] = help / Z->matrix[xx][xx];
		}
		iGS++;
		deleteMatrix(y);
		delete y;
		deleteMatrix(Z);
		delete Z;
		norm = calculateNorm(A, x, b);
		if (norm > lastNorm) {
			cout << "Metoda Gaussa-Seidela nie zbiega sie" << endl;
			file << "Metoda Gaussa-Seidela nie zbiega sie" << endl;
			return;
		}
	}

	deleteMatrix(D);
	delete D;
	deleteMatrix(L);
	delete L;
	deleteMatrix(U);
	delete U;

	double finalNorm = calculateNorm(A, x, b);
	long time = (clock() - start);

	cout << "Ilosc iteracji: " << iGS << endl;
	cout << "Norma: " << finalNorm << endl;
	cout << "Czas: " << time << " ms." << endl;

	file << "Ilosc iteracji: " << iGS << endl;
	file << "Norma: " << finalNorm << endl;
	file << "Czas: " << time << " ms." << endl;

}

void LU(MatrixStruct *A, MatrixStruct *x, MatrixStruct *b, fstream& file) {
	MatrixStruct *U = new MatrixStruct(A->m, A->n);
	MatrixStruct *L = new MatrixStruct(A->m, A->n);

	clock_t start = clock();

	fillMatrix(U, 0);
	fillMatrix(L, 1, 0, 0, 0);

	double help;
	for (int xx = 0; xx <= A->m - 1; xx++) {
		//U
		for (int yy = xx; yy <= A->n - 1; yy++) {
			help = 0;
			for (int zz = 0; zz <= xx - 1; zz++) {
				help += L->matrix[xx][zz] * U->matrix[zz][yy];
			}
			U->matrix[xx][yy] = A->matrix[xx][yy] - help;
		}
		//L
		for (int yy = xx + 1; yy <= A->m - 1; yy++) {
			help = 0;
			for (int zz = 0; zz <= xx - 1; zz++) {
				help += L->matrix[yy][zz] * U->matrix[zz][xx];
			}
			L->matrix[yy][xx] = (A->matrix[yy][xx] - help) / U->matrix[xx][xx];
		}
	}
	//Ax=b	A=LU
	//LUx=b
	//y=Ux
	//Ly=b
	//Ux=y
	cout << "LU:" << endl;
	file << "LU:" << endl;
	//Ly=b
	MatrixStruct *y = new MatrixStruct(L->m, 1);
	fillMatrix(y, 0);
	for (int xx = 0; xx < L->m; xx++) {
		help = b->matrix[xx][0];
		for (int zz = 0; zz < xx; zz++) {
			help -= (L->matrix[xx][zz]) * (y->matrix[zz][0]);
		}
		y->matrix[xx][0] = help / L->matrix[xx][xx];
	}
	deleteMatrix(L);
	delete L;
	//Ux=y
	for (int xx = 0; xx < U->m; xx++) {
		help = y->matrix[y->m - 1 - xx][0];
		for (int zz = 0; zz < xx; zz++) {
			help -= (U->matrix[y->m - 1 - xx][y->m - 1 - zz]) * (x->matrix[y->m - 1 - zz][0]);
		}
		x->matrix[y->m - 1 - xx][0] = help / U->matrix[y->m - 1 - xx][y->m - 1 - xx];
	}

	deleteMatrix(y);
	delete y;
	deleteMatrix(U);
	delete U;

	double finalNorm = calculateNorm(A, x, b);
	long time = (clock() - start);

	cout << "Norma: " << finalNorm << endl;
	cout << "Czas: " << time << " ms." << endl;

	file << "Norma: " << finalNorm << endl;
	file << "Czas: " << time << " ms." << endl;

}
int main() {
	//zad A
	int d = INDEKS % 10;
	int c = (INDEKS % 100 - d) / 10;
	int e = ((INDEKS % 1000) - 10 * c - d) / 100;
	int f = ((INDEKS % 10000) - 100 * e - 10 * c - d) / 1000;
	int N = 100 * 9 + 10 * c + d;
	//int N = 7;

	double a1 = 5 + e;
	double a2 = -1;
	double a3 = -1;

	MatrixStruct *A = new MatrixStruct(N, N);
	fillMatrix(A, a1, a2, a3, 0);

	MatrixStruct *b = new MatrixStruct(N, 1);
	for (int x = 0; x < b->m; x++) {
		for (int y = 0; y < b->n; y++) {
			b->matrix[x][y] = sin((x + 0)*(f + 1));
		}
	}

	MatrixStruct *x = new MatrixStruct(N, 1);
	fillMatrix(x, 0);

	fstream file("data.txt", ios::out);

	//zad B
	cout << "Zadanie B" << endl << endl;
	file << "Zadanie B" << endl << endl;

	Jacob(A, x, b, file);
	//x reset
	fillMatrix(x, 0);
	cout << endl;
	file << endl;

	GaussSeidel(A, x, b, file);
	//x reset
	fillMatrix(x, 0);
	cout << "-----------------------------------------" << endl;
	file << "-----------------------------------------" << endl;

	//zad C
	cout << "Zadanie C" << endl << endl;
	file << "Zadanie C" << endl << endl;
	a1 = 3;
	a2 = -1;
	a3 = -1;
	fillMatrix(A, a1, a2, a3, 0);

	Jacob(A, x, b, file);
	//x reset
	fillMatrix(x, 0);
	cout << endl;
	file << endl;

	GaussSeidel(A, x, b, file);
	//x reset
	fillMatrix(x, 0);
	cout << "-----------------------------------------" << endl;
	file << "-----------------------------------------" << endl;

	//zad D
	cout << "Zadanie D" << endl << endl;
	file << "Zadanie D" << endl << endl;

	LU(A, x, b, file);
	//x reset
	fillMatrix(x, 0);
	cout << "-----------------------------------------" << endl;
	file << "-----------------------------------------" << endl;

	deleteMatrix(A);
	delete A;
	deleteMatrix(x);
	delete x;
	deleteMatrix(b);
	delete b;

	//zad E
	a1 = 5 + e;
	a2 = -1;
	a3 = -1;
	cout << "Zadanie E" << endl << endl;
	file << "Zadanie E" << endl << endl;
	int unknown[] = { 100,500,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000 };
	//int unknown[] = { 100,500 };
	for (int i = 0; i < sizeof(unknown) / sizeof(unknown[0]); i++) {
		cout << "N = " << unknown[i] << endl;
		file << "N = " << unknown[i] << endl;

		MatrixStruct *A = new MatrixStruct(unknown[i], unknown[i]);
		fillMatrix(A, a1, a2, a3, 0);
		MatrixStruct *b = new MatrixStruct(unknown[i], 1);
		for (int x = 0; x < b->m; x++) {
			for (int y = 0; y < b->n; y++) {
				b->matrix[x][y] = sin((x + 0)*(f + 1));
			}
		}
		MatrixStruct *x = new MatrixStruct(unknown[i], 1);
		fillMatrix(x, 0);

		Jacob(A, x, b, file);
		//x reset
		fillMatrix(x, 0);
		cout << endl;
		file << endl;

		GaussSeidel(A, x, b, file);
		//x reset
		fillMatrix(x, 0);
		cout << endl;
		file << endl;

		LU(A, x, b, file);
		//x reset
		fillMatrix(x, 0);
		cout << endl;
		file << endl;

		deleteMatrix(A);
		delete A;
		deleteMatrix(x);
		delete x;
		deleteMatrix(b);
		delete b;
		cout << "*****************************************" << endl;
		file << "*****************************************" << endl;
	}


	file.close();

	return 0;
}
