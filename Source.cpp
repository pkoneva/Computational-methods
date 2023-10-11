#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void StraightWay(double** A, double **inv, double* b, int N) {
	double cos, sin, znam, rem1, rem2;
	int i, j, k;
	for (i = 0; i < N; i++) {
		for (j = i+1; j < N; j++) {
			znam = sqrt(A[i][i] * A[i][i] + A[j][i] * A[j][i]);
			cos = A[i][i] / znam;
			sin = A[j][i] / znam;
			
			for (k = i; k < N; k++) {
				rem1 = A[i][k];
			rem2 = A[j][k];
				A[i][k] = cos * rem1 + sin * rem2;
				A[j][k] = -sin * rem1 + cos * rem2;
				

			}
			for (k = 0; k < N; k++) {
				rem1 = inv[i][k];
				rem2 = inv[j][k];
				inv[i][k] = cos * rem1 + sin * rem2;
				inv[j][k] = -sin * rem1 + cos * rem2;
			}
			rem1 = b[i];
			rem2 = b[j];
			b[i]= cos * rem1 + sin * rem2;
			b[j]= -sin * rem1 + cos * rem2;
			
		}
	}
}


float chek(long double x) {
	float y = fabs(1 - x);
	if (y == 1)
		return 0;
	else return x;
}
int chek_b(double** A, double* b, int N) { //проверка столбца свободных коэф. на 0
	int i;
	int count0 = 0;

	for (i = N - 1; i >= 0; i--) {
		b[i] = float(b[i]);
		A[i][i] = float(A[i][i]);
		if (chek(b[i]) == 0 && chek(A[i][i]) == 0) count0++;
		if (chek(b[i]) != 0 && chek(A[i][i]) == 0) return -1;// у какого-то уравнения нет решения

	}
	if (count0 != 0) return 0;// есть уранвение вида 0=0,  в слу уравнений на 1 меньше, чем неизвестных, есть беск. мн-во решений
	return 1;
}

double Det(double** A, int N) {
	int i;
	double det = 1;
	for (i = 0; i < N; i++)
		det *= A[i][i];
	return det;
}

void GetX(double** M, int N, double* b,double* x) { // Обратный ход
	x[N - 1] = b[N - 1] / M[N - 1][N - 1];
	for (int i = N - 2; i >= 0; i--)
	{
		double rem = 0;
		for (int j = i + 1; j < N; j++)
		{
			rem = rem + M[i][j] * x[j]; //вычисляем, что нужно отнять от свободного коэффициента
		}
		x[i] = (b[i] - rem) / M[i][i];
	}
}


void Get_Inverse(double** M, double** inv, int N) {
	int i, j, k; 
	double coef;
	for (k = N - 1; k > 0; k--)
	{
		for (i = k - 1; i >= 0; i--)
		{
			coef = M[i][k] / M[k][k];
			for (j = N-1; j >=0; j--) {
				M[i][j] -= M[k][j] * coef;
				inv[i][j] -= inv[k][j] * coef;
			}
			
		}
		

	}
	for(k=0; k<N; k++)
		for (j = 0; j < N; j++) {
			inv[k][j] /= M[k][k];
		}

	printf("\n\nОбратная матрица:\n");
	for (j = 0; j < N; j++) {
		for (i = 0; i < N; i++)
			printf("%lf ", inv[j][i]);
		printf("\n");
	}

}



int main()
{

	int N;
	double** M;
	double* b;
	double** inverseM;
	int zero_b;
	double det;
	double* X;
	int i, j;

	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);
	scanf("%d", &N);

	M = (double**)calloc(N, sizeof( double*));
	inverseM = (double**)calloc(N, sizeof(double*));
	b = (double*)calloc(N, sizeof(double));

	for (i = 0; i < N; i++) {
		M[i] = (double*)calloc(N, sizeof(double));
		inverseM[i] = (double*)calloc(N, sizeof(double));
		for (j = 0; j < N; j++) {
			scanf("%Lf", &M[i][j]);
		}
		inverseM[i][i] = 1;
	}

	
	for (i = 0; i < N; i++)
		scanf("%Lf", &b[i]);

	StraightWay(M, inverseM, b, N);
	zero_b = chek_b(M, b, N);
	printf("\ndet=%Lf\n", det = chek(Det(M, N)));
	

	if (chek(det) == 0 && zero_b == -1) printf("Пустое множество решений. Обратной матрицы не существует.");
	else if (chek(det) == 0 && zero_b == 0)
		printf("Бесконечное множество решений. Обратной матрицы не существует. ");
	else {
		X = (double*)calloc(N, sizeof(double));
		GetX(M, N, b, X);
		printf("\nРешение:(");
		for (i = 0; i < N - 1; i++)
			printf("%Lf, ", X[i]);
		printf("%Lf)", X[N - 1]);
		free(X);
		Get_Inverse(M,inverseM,  N);
	}

	free(inverseM);
	free(M);
	free(b);
	return 0;
}