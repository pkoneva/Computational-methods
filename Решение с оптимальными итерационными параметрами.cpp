#include <fstream>
#include <iostream>
#include <math.h>
using namespace std;
const double PI = 3.1415926535897932384626433832795;
ifstream fin;
ofstream fout;
class Matrix {
	int n;//высота
	int m;//длина
	double** matr;
	double* vect;
public:
	Matrix() {
		this->n = 0;
		this->m = 0;
		this->matr = 0;
		this->vect = 0;
	}
	Matrix(int N) {
		if (N <= 0) throw "Неверный размер";
		this->n = N;
		this->m = N;
		this->matr = new double* [N];

		for (int i = 0; i < N; i++)
			this->matr[i] = new double[N];
	}
	Matrix(int N, int M) {
		if (N <= 0) throw "Неверный размер";
		if (M <= 0) throw "Неверный размер";
		this->n = N;
		this->m = M; if (M != 1) {
			this->matr = new double* [N];
			this->vect = 0;
			for (int i = 0; i < M; i++)
				this->matr[i] = new double[N];
		}
		else {
			this->matr = 0;
			this->vect = new double[N];
		}
	}
	~Matrix() {
		delete[] this->matr;
		delete[] this->vect;
	}
	Matrix(const Matrix& that) {
		this->n = that.n;
		this->m = that.m;
		if (that.m != 1) {
			this->matr = new double* [this->n];
			this->vect = 0;
			for (int i = 0; i < this->n; i++) {
				this->matr[i] = new double[this->m];
				for (int j = 0; j < this->n; j++) {
					this->matr[i][j] = that.matr[i][j];
				}
			}
		}
		else {
			this->matr = 0;
			this->vect = new double[this->n];
			for (int j = 0; j < this->n; j++) {
				this->vect[j] = that.vect[j];
			}
		}

	}
	Matrix operator+(const Matrix& that) const {
		if (this->n == that.n && this->m == that.m) {
			if (this->m != 1) {
				Matrix res(this->n, that.m);
				for (int i = 0; i < this->n; i++) {
					for (int j = 0; j < this->m; j++)
						res.matr[i][j] = this->matr[i][j] + that.matr[i][j];
				}
				return res;
			}
			else {
				Matrix res(this->n, that.m);
				for (int j = 0; j < this->n; j++)
					res.vect[j] = this->vect[j] + that.vect[j];
				return res;
			}
		}
		else throw "Неверный размер";
	}
	Matrix operator-(const Matrix& that) const {
		if (this->n == that.n && this->m == that.m) {
			if (this->m != 1) {
				Matrix res(this->n, that.m);
				for (int i = 0; i < this->n; i++) {
					for (int j = 0; j < this->m; j++)
						res.matr[i][j] = this->matr[i][j] - that.matr[i][j];
				}
				return res;
			}
			else {
				Matrix res(this->n, that.m);
				for (int j = 0; j < this->n; j++)
					res.vect[j] = this->vect[j] - that.vect[j];
				return res;
			}
		}
		else throw "Неверный размер";
	}
	Matrix operator*(const Matrix& that) const {
		if (this->m == that.n) {
			Matrix res(this->n, that.m);
			double sum = 0;
			if (that.m != 1) {
				for (int i = 0; i < this->n; i++) {
					for (int j = 0; j < that.m; j++) {
						for (int k = 0; k < this->m; k++)
							sum += this->matr[i][k] * that.matr[k][j];
						res.matr[i][j] = sum;
						sum = 0;
					}
				}
			}
			else {
				for (int j = 0; j < that.n; j++) {
					for (int k = 0; k < this->m; k++)
						sum += this->matr[j][k] * that.vect[k];
					res.vect[j] = sum;
					sum = 0;
				}
			}
			return res;
		}
		else throw "Неверный размер";

	}
	bool operator==(const Matrix& that) const {
		if (this->n == that.n && this->m == that.m) {
			if (that.m != 1) {
				for (int i = 0; i < this->n; i++) {
					for (int j = 0; j < that.m; j++) {
						if (this->matr[i][j] != that.matr[i][j]) return false;
					}
				}
				return true;
			}
			else {
				for (int j = 0; j < that.n; j++) {
					if (this->vect[j] != that.vect[j]) return false;
				}
				return true;
			}

		}
	}
	Matrix& operator=(const Matrix& that) {
		if (this != &that) {
			this->~Matrix();

			this->n = that.n;
			this->m = that.m;
			if (that.m != 1) {
				this->matr = new double* [this->n];
				this->vect = 0;
				for (int i = 0; i < this->n; i++) {
					this->matr[i] = new double[this->m];
					for (int j = 0; j < this->m; j++)
						this->matr[i][j] = that.matr[i][j];
				}
			}
			else {
				this->vect = new double[this->n];
				this->matr = 0;
				for (int j = 0; j < this->n; j++)
					this->vect[j] = that.vect[j];
			}
		}
		return *this;

	}
	void PutElem(int i, int j, double E) {
		if (i > this->n || j > this->m)throw "Неверный размер";
		if (this->m == 1) this->vect[i] = E;
		else this->matr[i][j] = E;
	}
	void SetElem() {
		if (this->m != 1) {
			for (int i = 0; i < this->n; i++)
				for (int j = 0; j < this->m; j++)
					fin >> this->matr[i][j];
		}
		else {
			for (int j = 0; j < this->n; j++)
				fin >> this->vect[j];
		}
	}
	void Transpose() {
		if (this->n != this->m) throw "Неверный размер";
		double rem;
		for (int i = 0; i < this->n; i++) {
			for (int j = i + 1; j < this->n; j++) {
				rem = this->matr[i][j];
				this->matr[i][j] = this->matr[j][i];
				this->matr[j][i] = rem;

			}
		}

	}
	void PrintMatr() {
		if (this->m != 1) {
			for (int i = 0; i < this->n; i++) {
				for (int j = 0; j < this->m; j++)
					fout << this->matr[i][j] << " ";
				fout << endl;
			}
		}
		else {
			for (int j = 0; j < this->n; j++)
				fout << this->vect[j] << " ";
		}
	}
	void Multiply_Scalar(double s) {
		if (this->m != 1) {
			for (int i = 0; i < this->n; i++) {
				for (int j = i + 1; j < this->n; j++) {
					this->matr[i][j] *= s;
				}
			}
		}
		else {
			for (int j = 0; j < this->n; j++) {
				this->vect[j] *= s;
			}
		}

	}
	double Norma_Vect_1() {
		double sum = 0;
		for (int i = 0; i < this->n; i++) {
			sum += fabs(this->vect[i]);
		}
		return sum;
	}

};


int main() {
	int N;
	double L1, L2;
	double t;
	double k, p;
	double e = 0.000001;
	fin.open("input.txt");
	fout.open("output.txt");
	fin >> N;
	fin >> L1 >> L2;
	Matrix X(N, 1);
	X.PutElem(0, 0, 1);
	for (int i = 1; i < N; i++) {
		X.PutElem(i, 0, 0);
	}
	Matrix A(N);
	A.SetElem();
	Matrix B(N, 1);
	B.SetElem();

	Matrix Rem(N, 1);
	p = (1 - sqrt(L1 / L2)) / (1 + sqrt(L1 / L2));//вычисление длины цикла итераций
	k = log(2 / e) / log(1 / p);
	k = abs(k) + 1;
	fout << k;

	int nk = abs(log(k) / log(2));//новый размер массива для работы с оптимальной последовательностью
	k = nk;
	nk = 1;
	for (int i = 0; i <= k; i++) {
		nk = 2 * nk;
	}
	fout << endl << nk << endl;


	int* NK = new int[nk];//массив номеров параметров
	for (int i = 0; i < nk; i++) NK[i] = i;
	int* a = new int[nk];
	int m = 0; int s = 1;
	while (s < nk / 2) {
		int y = 0;
		int l = 0, r = nk - s;
		while (l < r - s+1) {
			for (int i = 0; i < s; i++) a[y++] = NK[l + i];
			for (int i = 0; i < s; i++) a[y++] = NK[r + i];
			l += s;
			r -= s;
		}
		m++;
		s *= 2;
		for (int i = 0; i < nk; i++) {
			NK[i] = a[i];
		}
	}
	delete(a);

	double* T = new double[nk + 1];//массив параметров
	for (int i = 1; i < nk + 1; i++) {
		T[i] = 2 / ((L2 + L1) + (L2 - L1) * fabs(cos((2 * i - 1) * PI / (2 * k))));
	}
	Matrix help(N, 1);
	help = A * X - B;
	while (help.Norma_Vect_1() > e) {
		int count = 0;
		for (int i = 1; i <= k; i++) {
			help.Multiply_Scalar(T[NK[i]]);
			Rem = X - help;
			X = Rem;
			help = A * X - B;
		}
	}
	fout << "\n";
	X.PrintMatr();
	fin.close();
	fout.close();
	delete(NK);
	delete(T);
	return 0;
}