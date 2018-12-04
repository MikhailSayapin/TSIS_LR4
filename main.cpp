#include <iostream>
#include <cmath>

using namespace std;

double func(double x) {
	/*return -sqrt(x)*sin(x)*sin(5*x);*/
	return (cos(x) + log10(x))*sin(5 * x);
}

double random(double x1,double x2) {
	int num = 100000000;
	return fabs((double)((rand()*num) % (int)((x2 - x1)*num) + 1) / num) + x1;
}


int main() {
	double a = 7;
	double b = 10;
	double p;
	double temp = 10000000;
	double x1 = random(a,b);
	double x2;

	while (temp > 0.05) {
		x2 = random(a, b);
		if (func(x2) < func(x1)) {
			x1 = x2;
		}
		else if (func(x2) >= func(x1)) {
			p = rand() % 100;
			p /= 100;
			if (p < exp(-(func(x2) - func(x1) / temp))) {
				x1 = x2;
			}
			
		}
		temp *= 0.95; //0.95
		cout << temp << endl;
	}

	cout << func(x1) << endl;
	cout << func(x2) << endl;



	return 0;
}
