#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

const int N = 31;
const double h = 1. / (N - 1);
const double PI = 3.14159265358979;

const double tay = 0.1; // шаг по времени тау
const double e = 0.01; //заданная точность в относительных единицах
const double alfa = 1.; // заданная альфа для давления

ofstream a("fout.txt");

void initialPre(int n, double A[]) {
	for (int j = 0; j < n; j++) {
		A[j] = alfa * j * h * (1 - j * h);
	}
}

void null(int n, double A[]) {
	for (int j = 0; j < n; j++) {
		A[j] = 0.;
	}
}

void output(int n, double A[]) {
	for (int j = 0; j < n; j++) {
		cout<<floor(A[j] * 10000) / 10000<<" ";
	}
	cout<<endl;
}

void outputFile(int n, double A[], ofstream& a) {
	for (int j = 0; j < n; j++) {
		a<<floor(A[j] * 10000) / 10000<<" "<<endl;
	}
	a<<endl;
}

void lamda(double n, double A[]) {
	for (int j = 0; j < n; j++) {
		A[j] = 4. / (h * h) * sin(j * PI * h / (2.)) * sin(j * PI * h / (2.));
	}
}

void furie(int n, double A[], double S[]) {
	double sum, norm = sqrt(2. / (n - 1));
	for (int i = 0; i < n; i++) {
		sum = 0;
		for (int j = 0; j < n; j++) {
			sum += A[j] * sin(fmod(PI * i * j / (n - 1), 2 * PI)); //fmod вычисляет остаток от деления, результат с плав. точкой
		}
		S[i] = sum * norm;
	}
}

void fcsp1(int n, double CS[], double CSP1[], double pre[], double l[], double t) {
	for (int i = 0; i < n; i++) {
		CSP1[i] = (pre[i] * t + CS[i]) / (1 + l[i] * t);
	}
}

void assignm(int n, double A[], double B[]) {
	for (int i = 0; i < n; i++) {
		B[i] = A[i];
	}
}

void outputmas(int n, double cs[], double ws[], double csp1[], double wsp1[]) {
	cout<<"Cs: ";
	output(N, cs);
	cout<<"Ws: ";
	output(N, ws);
	cout<<"Cs+1: ";
	output(N, csp1);
	cout<<"Ws+1: ";
	output(N, wsp1);
	cout<<endl;
}

void accuracyy(int n, double& e1, double A[], double B[], int& k) {
	double e2 = 0; //точность в точке i
	double e2p1 = 0; //точность в точке i+1
	for (int i = 1; i < n-1; i++) {
		e2p1 = 2 * ( abs(B[i]- A[i]) ) / ( abs(A[i]) + abs(B[i]) );
		if (e2p1 > e2) {
			e2 = e2p1;
			k = i;
		}
	}
	e1 = e2;
}

void initialhm(int n, double hm[]) {
	for (int i = 0; i < n; i++) {
		hm[i] = i * h;
	}
}

void analiticResh(int n, double resh[]) {
	double x = 0;
	for (int i = 0; i < n; i++) {
		x = i * h;
		resh[i] = (-1) * (-x / 12 + x * x * x / 6 - x * x * x * x / 12);
	}
}

void statProgib() {

	double cs[N], ws[N], csp1[N], wsp1[N], pre[N], pref[N], l[N], hm[N];

	initialPre(N, pre); //инициализация массива давления
	furie(N, pre, pref); // нахождение собственных значений давления
	initialhm(N, hm); //инициализация массива шагов
	lamda(N, l); //инициализация лямбд

	null(N, cs); //обнуление массива Cs
	null(N, ws); //обнуление массива Ws
	null(N, csp1); //обнуление массива Cs+1
	null(N, wsp1); //обнуление массива Ws+1

	cout<<"p: ";
	output(N, pre); // вывод массива давления
	cout<<"pf: ";
	output(N, pref); // вывод массива давления
	cout<<"l: ";
	output(N, l); //вывод массива лямбд
	cout<<endl;

	outputmas(N, cs, ws, csp1, wsp1);

	int s = 0; //номер текущего шага
	double e1 = 1; // точность текущего шага
	int k = 0; //номер элемента по которому считает точность

	while (e1 > e) {
		assignm(N, wsp1, ws); // Ws=Ws+1
		furie(N, ws, cs); // находим Cs преобр. Фурье Ws
		fcsp1(N, cs, csp1, pref, l, tay); // находим Cs+1 из Cs
		furie(N, csp1, wsp1); // находим Ws+1 преобр. Фурье Cs+1

		outputmas(N, cs, ws, csp1, wsp1); //выводим все массивы

		accuracyy(N, e1, ws, wsp1, k); //вычисляем точность в относительных величинах

		cout<<"e1= "<<e1<<" k="<<k<<" s="<<s<<endl<<endl; //вывод информации в консоль
		s += 1;
	}

	a<<"pre:"<<endl; //блок выводов в фаил
	outputFile(N, pre, a);
	a<<"wsp1:"<<endl;
	outputFile(N, wsp1, a);
	a<<"hm:"<<endl;
	outputFile(N, hm, a);

	double resh[N]; //нахождение и вывод в фаил аналитического решения
	analiticResh(N, resh);
	a<<"resh:"<<endl;
	outputFile(N, resh, a);

}

int main() {

	statProgib();
	a.close();
	system("pause");

}

