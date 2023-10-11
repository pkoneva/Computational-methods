#include <fstream>
#include <iostream>
#include <math.h>
#include <algorithm>
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

	Matrix& Transpose() {
		double rem;
		if (this->n == this->m) {

			for (int i = 0; i < this->n; i++) {
				for (int j = i + 1; j < this->n; j++) {
					rem = this->matr[i][j];
					this->matr[i][j] = this->matr[j][i];
					this->matr[j][i] = rem;

				}
			}
		}
		Matrix T(m, n);
		for (int i = 0; i < n; i++) {
			for (int j = i + 1; j < m; j++) {
				rem = this->matr[i][j];
				T.matr[j][i] = rem;

			}
		}
		return T;
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
	double Scal_proizv(const Matrix& y) {
		if (n != y.n || m != 1 || y.m != 1) throw "Неверные значения";
		double res = 0;
		for (int i = 0; i < n; i++)
			res += vect[i] * y.vect[i];
		return res;

	}
	void Normalize_vekt() {
		if (m != 1) throw "Неверные значения";
		double norma = Norma_Vect_1();
		if (norma != 0) {
			for (int i = 0; i < n; i++)
				vect[i] /= norma;
		}
		return;
	}

	void Metod_Givensa() {
		//double cos, sin, znam, rem1, rem2;
		//int i, j, k;
		//for (i = 0; i < n - 2; i++) {
		//	for (j = i + 2; j < n; j++) {
		//		znam = sqrt(matr[i + 1][i] * matr[i + 1][i] + matr[j][i] * matr[j][i]);
		//		cos = matr[i + 1][i] / znam;
		//		sin = matr[j][i] / znam;
		//		matr[i][j] = 0/*cos*/;
		//		matr[j][i] = 0/*sin*/;
		//		matr[i][i + 1] = znam;
		//		matr[i + 1][i] = znam;

		//		for (k = i + 1; k < n; k++) {
		//			//rem1 = matr[k][i+1];
		//			//rem2 = matr[k][j];
		//			//matr[k][i+1] = cos * rem1 + sin * rem2;
		//			//matr[k][j] = -sin * rem1 + cos * rem2;
		//			///*rem1 = matr[i+1][k];
		//			//rem2 = matr[j][k];
		//			//matr[i+1][k] = cos * rem1 + sin * rem2;
		//			//matr[j][k] = -sin * rem1 + cos * rem2;*/
		//			//matr[i + 1][k] = matr[k][i + 1];
		//			//matr[j][k] = matr[k][j];


		//			rem1 = matr[k][i + 1];
		//			rem2 = matr[k][j];
		//			matr[k][i + 1] = cos * rem1 + sin * rem2;
		//			matr[k][j] = -sin * rem1 + cos * rem2;
		//			matr[i + 1][k] = matr[k][i + 1];
		//			matr[j][k] = matr[k][j];
		//		}

		//	}
		//}
Matrix res(n);
		for (unsigned int i = 0; i < n /*- 1*/; i++)
		{
			for (unsigned int j = i + 2; j < n; j++)
			{
				if (matr[i][i] - matr[j][j] == 0) {
					fout << "матрица не может быть приведена к трёхдиагональному виду" << endl;
					return;
				}
				double t = 2 * matr[i][j] / (matr[i][i] - matr[j][j]);
				double phi = 0.5 * atan(t);
				double c = cos(phi);
				double s = sin(phi);

				/*float bii = c * c * matr[i][i] + 2 * c * s * matr[i][j] + s * s * matr[j][j];
				float bij = s * c * (matr[j][j] - matr[i][i]) + matr[i][j] * (c * c - s * s);
				float bjj = s * s * matr[i][i] + c * c * matr[j][j] - 2 * c * s * matr[i][j];
				float bji = bij;

				matr[i][i] = bii;
				matr[i][j] = bij;
				matr[j][i] = bji;
				matr[j][j] = bjj;*/
				/*double c, s, znam, rem1, rem2;
				znam = sqrt(matr[i][i] * matr[i][i] + matr[j][i] * matr[j][i]);
					c = matr[i][i] / znam;
				s = matr[j][i] / znam;*/
				Matrix B(n);
				for (int k = 0; k < n; k++) {
					for (int l = 0; l < n; l++) B.matr[k][l] = 0;
					B.matr[k][k] = 1;
				}
				B.matr[i][i] = c;
				B.matr[j][j] = c;
				B.matr[i][j] = s;
				B.matr[j][i] = -s;
				B.PrintMatr();
				res = B * *this;
				res.PrintMatr();
				B.matr[i][j] = -s;
				B.matr[j][i] = s;
				B.PrintMatr();
				*this = res;
				res =  *this * B;
				res.PrintMatr();
				*this = res;
				/*matr[j][i] = 0;
				matr[i][j] = 0;*/
			}
			//for (unsigned int i = 0; i < n /*- 1*/; i++)
			//{
			//	for (unsigned int j = i + 2; j < n; j++) {
			//		matr[j][i] = 0;
			//		matr[i][j] = 0;
			//	}
			//}
		}

		/*for (int i = 0; i < n - 1; i++)
		{
			for (int j = i + 1; j < n; j++)
			{
				if (matr[i][i] == 0 && matr[j][i] == 0)

					int y;
				else

				{
					double c = matr[i][i] / sqrt(matr[i][i] * matr[i][i] + matr[j][i] * matr[j][i]);

					double s = matr[j][i] / sqrt(matr[i][i] * matr[i][i] + matr[j][i] * matr[j][i]);

					for (int m = i; m < n; m++)
					{
						double h = matr[i][m], g = matr[j][m];

						matr[i][m] = h * c + g * s;

						matr[j][m] = g * c - h * s;

					}
				}
			}
		}*/
	}

	void Metod_xays() {
		int l, k, j, i;
		float scale, hh, h, g, f;
		Matrix e(n, 1), d(n, 1);
		e.vect[0] = 0;
		/* Проход по стадиям процесса редукции */
		for (i = n-1; i >= 1; i--) {
			l = i - 1; h = scale = 0.;
			/* сложный процесс везде, кроме последней стадии */
			if (l > 1) {
				/* вычислить шкалу */
				for (k = 0; k < l; k++) scale += fabs(matr[i][k]);
				/* малая величина шкалы -> пропустить преобразование */
				if (scale == 0.) e.vect[i] = matr[i][l];
				else {
					/* отмасштабировать строку и вычислить s2 в h */
					for (k = 0; k < l; k++) {
						matr[i][k] /= scale; h += matr[i][k] * matr[i][k];
					}
					/* вычислить вектор u */
					f = matr[i][l];
					g = (f >= 0. ? -sqrt(h) : sqrt(h));
					e.vect[i] = scale * g; h -= f * g;
					/* записать u на место i-го ряда a */
					matr[i][l] = f - g;
					/* вычисление u/h, Au, p, K */
					f = 0.;
					for (j = 0; j < l; j++) {
						/* следующая инструкция не нужна, если не требуются вектора,
							  она содержит загрузку u/h в столбец a */
			/*			matr[j][i] = matr[i][j] / h;*/
						/* сформировать элемент Au (в g) */
						g = 0.;
						for (k = 0; k < j; k++) g += matr[j][k] * matr[i][k];
						for (k = j; k < l; k++) g += matr[k][j] * matr[i][k];
						/* загрузить элемент p во временно неиспользуемую область e */
						e.vect[j] = g / h;
						/* подготовка к формированию K */
						f += e.vect[j] * matr[i][j];
					}
					/* Сформировать K */
					hh = f / (h + h);
					for (j = 0; j < l; j++) {
						/* Сформировать q и поместить на место p (в e) */
						f = matr[i][j]; e.vect[j] = g = e.vect[j] - hh * f;
						/* Трансформировать матрицу a */
						for (k = 0; k < j; k++) matr[j][k] -= (f * e.vect[k] + g * matr[i][k]);
					}
				}
			}
			else e.vect[i] = matr[i][l];
			d.vect[i] = h;
		}
		/* если не нужны собственные вектора, опустите следующую инструкцию */
	/*	d.vect[0] = 0.;*/
		/* эту опускать не надо */
		e.vect[0] = 0.;
		/* Все содержание цикла, кроме одной инструкции, можно опустить, если не
		   требуются собственные вектора */
		for (i = 0; i < n; i++) {
			//l = i - 1;
			///* этот блок будет пропущен при i=1 */
			//if (d.vect[i] != 0.) {
			//	for (j = 0; j < l; j++) {
			//		g = 0.;
			//		/* формируем PQ, используя u и u/H */
			//		for (k = 0; k < l; k++) g += matr[i][k] * matr[k][j];
			//		for (k = 0; k < l; k++) matr[k][j] -= g * matr[k][i];
			//	}
			//}
			/* эта инструкция остается */
			d.vect[i] = matr[i][i];
			/* ряд и колонка матрицы a преобразуются к единичной, для след. итерации */
			/*matr[i][i] = 0.;
			for (j = 0; j < l; j++) matr[j][i] = matr[i][j] = 0.;*/
		}

		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				matr[i][j] = 0;
				if (i == j) {
					matr[i][i] = d.vect[i];
					if (i != n - 1) {
						/*matr[i - 1][i] = e.vect[i + 1];*/
						matr[i][i + 1] = e.vect[i + 1];
						matr[i + 1][i] = matr[i][i - 1];
					}
					/*else {
						matr[i - 1][i] = matr[i - 1][i];
						matr[i][i - 1] = e.vect[i + 1];
					}*/
				}
				
			}
		}
	}

	double Norma_matr_1() {
		double sum = 0, max = 0;
		for (int j = 0; j < n; j++) {
			for (int i = 0; i < n; i++) {
				sum += fabs(matr[i][j]);
			}
			max = (max > sum ? max : sum);
			sum = 0;
		}
		return max;
	}
	double y1() {
		double y1;
		double help1, help2, help3 = 0;
		help1 = matr[0][0] - fabs(matr[0][1]);
		help2 = matr[n - 1][n - 1] - fabs(matr[n - 1][n - 2]);
		for (int i = 1; i < n - 1; i++) {
			double minrem = matr[i][i] - fabs(matr[i][i + 1]) - fabs(matr[i][i - 1]);
			help3 = (help3 < minrem ? help3 : minrem);
		}y1 = min(help3, min(help1, help2));
		return y1;
	}
	double y2() {
		double y2;
		double help1, help2, help3 = 0;
		help1 = matr[0][0] + fabs(matr[0][1]);
		help2 = matr[n - 1][n - 1] + fabs(matr[n - 1][n - 2]);
		for (int i = 1; i < n - 1; i++) {
			double maxrem = matr[i][i] + fabs(matr[i][i + 1]) + fabs(matr[i][i - 1]);
			help3 = (help3 > maxrem ? help3 : maxrem);
		}
		y2 = max(help3, max(help1, help2));
		return y2;
	}
	void bisection_method(double y1, double y2, double a, double b, int i) {
		static int count = 1;
			/*const int countrem = count;*/
			//help1 = matr[0][0] - fabs(matr[0][1]);
			//help2= matr[n-1][n-1] - fabs(matr[n-1][n-2]);
			//for (int i = 1; i < n - 1; i++) {
			//	double minrem= matr[i][i] - fabs(matr[i][i+1]) - fabs(matr[i+1][i+2]);
			//	help3 = (help3 > minrem ? help3 : minrem);
			//}
			//y1 = min(help3, min(help1, help2));
			//help1 = matr[0][0] + fabs(matr[0][1]);
			//help2 = matr[n - 1][n - 1] + fabs(matr[n - 1][n - 2]);
			//for (int i = 1; i < n - 1; i++) {
			//	double maxrem = matr[i][i] + fabs(matr[i][i + 1]) + fabs(matr[i + 1][i + 2]);
			//	help3 = (help3 > maxrem ? help3 : maxrem);
			//}
			//y2 = max(help3, max(help1, help2));*/
			//if (count == 0) {
			//	y1 = this->y1();
			//	y2 = this->y2();
			//}
		double l;
		/*const double rem1 = y1, rem2 = y2;*/
		double e = 0.00001;
		Matrix Anext(n);
		while (count <= n) {
			/*double reml;
			l = (rem1 + rem2) / 2;
			reml = l;*/
			if (/*countrem*/i < count) {
				y1 = a;
				y2 = b;
				i = count;
			}
			l = (y1 + y2) / 2;
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					if (i == j) Anext.matr[i][j] = matr[i][j] - l;
					else Anext.matr[i][j] = matr[i][j];
				}
			}
			int otr = 0, pol = 0;
			double dnext, d = matr[0][0] - l;
			for (int i = 1; i < n; i++) {
				if (d < 0) otr++;
				else if (d > 0) pol++;
				else;
				dnext = (matr[i][i] - l) - (matr[i][i - 1] * matr[i][i - 1] / d);
				d = dnext;
			}
			if (d < 0) otr++;
			else if (d > 0) pol++;
			else;
			if (y2 - y1 < e) {
				fout << (y1+y2)/2 << endl;
				count++;
				return;
			}
			/*if (otr >= count) {
				rem2 = l; continue;

			}*/


			if (otr >= count) {
				this->bisection_method(y1, l,a, b,  i);
				if (otr < count) break;

			}
			/*else {
				Anext.bisection_method(l, y2);
				return;
			}*/
			else/* (otr < count)*/ {
				this->bisection_method(l, y2,a, b, i);
				if (otr >= count) break;
			}
			/*
						if (otr >= count) {

						}*/
						/*	if (otr != 0) {double on = y2 - y1;
											if (on < e) fout << on << endl;
										while (otr != 1) {
											Anext.bisection_method(y1, l);
										}
										if (otr == 1) {

											while (on > e) {
												Anext.bisection_method(y1, l);

											}
										}
									}if (pol != 0) {double on = y2 - y1;
											if (on < e) fout << on << endl;
										while (pol != 1) {
											Anext.bisection_method(l, y2);
										}
										if (pol == 1) {

											while (y2 - y1 > e) Anext.bisection_method(l, y2);
										}
									}*/
		}
	}
};




int main() {
	int N;
	double own_num;
	fin.open("input.txt");
	fout.open("output.txt");
	fin >> N;
	/*Matrix X(N, 1);
	Matrix Xnext(N, 1);
	for (int i = 0; i < N; i++) {
		X.PutElem(i, 0, i + 1);
	}*/
	Matrix A(N);
	A.SetElem();
	A.Metod_Givensa();
	const double y1 = A.y1();
	const double y2 = A.y2();
	A.bisection_method(y1, y2, y1, y2, 1);
	A.PrintMatr();
	fin.close();
	fout.close();
	return 0;
}

