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
		this->m = M;
		if (M != 1) {
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
		if (n == 1) return;
		if(m!=1) delete[] this->matr;
		if(m==1) delete[] this->vect;
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
			if (that.m != 1 && m!=1) {
				for (int i = 0; i < this->n; i++) {
					for (int j = 0; j < that.m; j++) {
						for (int k = 0; k < this->m; k++)
							sum += this->matr[i][k] * that.matr[k][j];
						res.matr[i][j] = sum;
						sum = 0;
					}
				}
			}
			else if (that.m == 1 && m!=1){
				for (int j = 0; j < that.n; j++) {
					for (int k = 0; k < this->m; k++)
						sum += this->matr[j][k] * that.vect[k];
					res.vect[j] = sum;
					sum = 0;
				}
			}
			else {
				for (int j = 0; j < n; j++) {
					
					for (int k = 0; k < that.m; k++) {
						sum += that.matr[0][k] * vect[j];
						res.matr[j][k] = sum;
						
						sum = 0;
					}
					
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
	double GetElem(int i, int j) {
		if (i > this->n || j > this->m)throw "Неверный размер";
		return matr[i][j];
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
	void Transpose_n_m(Matrix & Y) {
		if (m != 1) {
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < m; j++) {
					Y.matr[j][i] = matr[i][j];
				}
			}
		}
		if (m == 1) {
			for (int i = 0; i < n; i++) {
					Y.matr[0][i] = vect[i];
			}
		}
		else throw "Неверный размер";
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
	double Norma_Vect_2() {
		double sum = 0;
		for (int i = 0; i < this->n; i++) {
			sum += this->vect[i]* this->vect[i];
		}
		sum = sqrt(sum);
		return sum;
	}
	double Scal_proizv(const Matrix& y) {
		if (n != y.n || m != 1 || y.m != 1) throw "Неверные значения";
		double res = 0;
		for (int i = 0; i < n; i++)
			res += vect[i] * y.vect[i];
		return res;

	}
	void Normalize_vekt() {
		if (m != 1) throw "Неверные значения";
		double norma = Norma_Vect_2();
		if (norma != 0) {
			for (int i = 0; i < n; i++)
				vect[i] /= norma;
		}
		return;
	}
	Matrix& StraightWay_Metod_Vracsheniy(Matrix& inv) {
		long double cos, sin, znam, rem1, rem2;
		int i, j, k;
		for (i = 0; i < n; i++) {
			for (j = i + 1; j < n; j++) {
				znam = sqrt(matr[i][i] * matr[i][i] + matr[j][i] * matr[j][i]);
				cos = matr[i][i] / znam;
				sin = matr[j][i] / znam;

				for (k = i; k < n; k++) {
					rem1 = matr[i][k];
					rem2 = matr[j][k];
					matr[i][k] = cos * rem1 + sin * rem2;
					matr[j][k] = -sin * rem1 + cos * rem2;


				}
				for (k = 0; k < n; k++) {
					rem1 = inv.matr[i][k];
					rem2 = inv.matr[j][k];
					inv.matr[i][k] = cos * rem1 + sin * rem2;
					inv.matr[j][k] = -sin * rem1 + cos * rem2;
				}
			}
		}
		return inv;
	}
	Matrix& Get_Inverse(Matrix& inv) {
		inv = StraightWay_Metod_Vracsheniy(inv);
		int i, j, k;
		long double coef;
		double e = 0.0000001;
		for (k = n - 1; k > 0; k--)
		{
			for (i = k - 1; i >= 0; i--)
			{
				coef = matr[i][k] / matr[k][k];
				for (j = n - 1; j >= 0; j--) {
					matr[i][j] -= matr[k][j] * coef;
					inv.matr[i][j] -= inv.matr[k][j] * coef;
				}

			}


		}
		for (k = 0; k < n; k++)
			for (j = 0; j < n; j++) {
				inv.matr[k][j] /= matr[k][k];
			}

		return inv;
	}
	double Find_own_number(Matrix& X) {
		static int i = n;
		Matrix Xnext(n, 1);
		long double lim_prev = 1;
		long double lim_next = 1;
		
		int count = 0;
		double e = 0.0000001;
		do {
			Xnext = X;
			Xnext.Normalize_vekt();
			Xnext = *this * Xnext;
			lim_prev = X.Norma_Vect_2();
			lim_next = Xnext.Norma_Vect_2();
			X = Xnext;
		} while (fabs(lim_next - lim_prev) > e);
		X.Normalize_vekt();
		fout << "Собственный вектор для" << i << "-ого собственного числа:" << endl;
		Xnext.PrintMatr();
		fout << endl << "Собственное число " << i << ":" << lim_next << endl;
		i--;
		return lim_next;
	}

};


int main() {
	int N;
	double e = 0.0000001;
	double own_num;
	fin.open("input.txt");
	fout.open("output.txt");
	fin >> N;
	Matrix X(N, 1);
	Matrix Xnext(N, 1);
	for (int i = 0; i < N; i++) {
		X.PutElem(i, 0, i);
	}
	Matrix A(N);
	A.SetElem();
	own_num = A.Find_own_number(X);
	Matrix B(N);
	Matrix Anext(N);
	Matrix Y(1, N);
	Matrix Q(N);
	for (int k = N - 1; k > 0; k--) {
		
		X.Transpose_n_m(Y);
		Q= X * Y;
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				if (i == j) B.PutElem(i, i, 1 - Q.GetElem(i, j));
				else B.PutElem(i, j, -Q.GetElem(i, j));
			}
		}
		Anext = B * A * B;
		for (int i = 0; i < N; i++) {
			X.PutElem(i, 0, i + 1);
		}
		Anext.Find_own_number(X);
		A = Anext;
	}


	fin.close();
	fout.close();
	return 0;
}
