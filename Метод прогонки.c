#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


void Straight_way(int N, double* A, double* B, double* C, double* d, double* Ac, double* Bc, double* Cc) {
	int i;
	Cc[0] = B[0];
	Bc[0] = d[0] / Cc[0];
	Ac[0] = -C[0] / Cc[0];
	for (i = 1; i < N-1; i++) {
		Cc[i] = B[i] + A[i] * Ac[i - 1];
		Bc[i] = (d[i] - A[i] * Bc[i - 1]) / Cc[i];
		Ac[i] = -C[i] / Cc[i];
	}
	Cc[N - 1] = B[N - 1] + A[N - 1] * Ac[N - 2];
	Bc[N-1] = (d[N-1] - A[N-1] * Bc[N - 2]) / Cc[N-1];
}

void Opposit_way(int N, double *X, double* Ac, double* Bc) {
	int i;
	X[N - 1] = Bc[N - 1];
	for (i = N - 2; i >= 0; i--) {
		X[i] = Ac[i] * X[i + 1] + Bc[i];
	}
}

int Chek_Matrix_Predominance_Condition(int N, double* A, double* B, double* C) {
	if(fabs(B[0]) < fabs(C[0])) {
		printf("Решение невозможно получить методом прогонки. Не выполняется условие диагонального преобладания.");
		return -1;
	}
	for (int i = 1; i < N-1; i++) {
		if (fabs(B[i]) < fabs(C[i]) + fabs(A[i]) || fabs(B[i]) <= fabs(A[i])) {
			printf("Решение невозможно получить методом прогонки. Не выполняется условие диагонального преобладания.");
			return -1;
		}
	}
	if (fabs(B[N-1]) < fabs(A[N-1])) {
		printf("Решение невозможно получить методом прогонки. Не выполняется условие диагонального преобладания.");
		return -1;
	}
	return 1;
}

 double Det(int N, double* A, double* B, double* C) {
	 if (N == 1) return B[0];
	 if (N == 0) return 1;
	 return B[N - 1] * Det(N - 1, A, B, C) - C[N - 2] * A[N - 1] * Det(N - 2, A, B, C);
}
 int chek_d(double* d, int N, double* A, double* B, double* C) {
	 int count0 = 0;
	 for (int i = N - 1; i >= 0; i--) {
		 if (d[i] != 0 && A[i] ==0&& B[i] == 0 && C[i] == 0) return -1;
		 if (d[i] == 0 && A[i]==0 && B[i] == 0 && C[i] == 0) count0++;
	 }
	 if (count0 != 0) return 0;
	 return 1;
 }


 void Inverse(int N, double* A, double* B, double *C){
	 int i, j;
	 double **inv= (double**)calloc(N, sizeof(double*));
	 double *e= (double*)calloc(N, sizeof(double));
	 double *Ac = (double*)calloc(N, sizeof(double));
	 double* Bc = (double*)calloc(N, sizeof(double));
	 double* Cc = (double*)calloc(N, sizeof(double));
	 for (i = 0; i < N; i++) 
		 inv[i] = (long double*)calloc(N, sizeof(long double));
	 for (i = 0; i < N; i++) {
		 e[i] = 1;
		 Straight_way(N, A, B, C, e, Ac, Bc, Cc);
		 Opposit_way(N, inv[i], Ac, Bc);
		 e[i] = 0;
	 }
	 for (j = 0; j < N; j++) {
		 for (i = 0; i < N; i++)
			 printf("%lf ", inv[i][j]);
		 printf("\n");
	 }
	 free(e);
	 free(inv);
	 free(Ac);
	 free(Bc);
	 free(Cc);
 }


int main()
{

	int N;
	double* A,* B, *C;
	double* Acoef, *Bcoef, *Ccoef;
	double* d;
	int true, ch_d;
	double det;
	double* X;
	int i;
	
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);
	scanf("%d", &N);

	A = (double*)calloc(N, sizeof(double)); // Нижняя диаг.
	B = (double*)calloc(N, sizeof(double)); // Главная
	C = (double*)calloc(N, sizeof(double)); // Верняя
	Acoef = (double*)calloc(N, sizeof(double));
	Bcoef = (double*)calloc(N, sizeof(double));
	Ccoef = (double*)calloc(N, sizeof(double));
	d = (double*)calloc(N, sizeof(double));

	for (i = 1; i < N; i++) 
		scanf("%lf", &A[i]);
	for (i = 0; i < N; i++)
		scanf("%lf", &B[i]);
	for (i = 0; i < N-1; i++)
		scanf("%lf", &C[i]);
	for (i = 0; i < N; i++)
		scanf("%Lf", &d[i]);

	ch_d = chek_d(d, N, A, B, C);
	printf("\ndet=%Lf\n", det = Det(N, A, B, C));
	if (det == 0 && ch_d == -1) printf("Пустое множество решений. Обратной матрицы не существует.");
	else if (det == 0 && ch_d == 0)
		printf("Бесконечное множество решений. Обратной матрицы не существует.");
	else {

		true = Chek_Matrix_Predominance_Condition(N, A, B, C);
		if (true == -1) return 0;

		Straight_way(N, A, B, C, d, Acoef, Bcoef, Ccoef);
		X = (double*)calloc(N, sizeof(double));
		Opposit_way(N, X, Acoef, Bcoef);
		printf("\nРешение:(");
		for (i = 0; i < N - 1; i++)
			printf("%.6lf, ", X[i]);
		printf("%.6lf)", X[N - 1]);
		free(X);
		printf("\n\nОбратная матрица:");
		Inverse(N, A, B, C);
	}

	free(A);
	free(B);
	free(C);
	free(Acoef);
	free(Bcoef);
	free(Ccoef);
	free(d);
	return 0;
}
