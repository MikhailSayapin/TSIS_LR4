#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

const int xmin = 0;
const double xmax = 3.14159;
const int K = 100;
const int M = 1;
const int L = 10 + 1;
const double P = 0.95;
const double eps = 0.01;
const int N = (log(1 - P) / log(1 - (eps / (xmax - xmin))));

double random(double x1, double x2)
{
	int num = 100000000;
	return fabs((double)((rand()*num) % (int)((x2 - x1)*num) + 1) / num) + x1;
}

int main()
{
	double xk[K + 1];

	for (int k = 0; k <= K; k++)
	{
		xk[k] = xmin + (k * (xmax - xmin)) / 100;
	}

	double f_shum[K + 1];
	for (int k = 0; k <= K; k++)
	{
		f_shum[k] = sin(xk[k]) + 0.5;
		f_shum[k] += random(-0.25, 0.25);
	}

	double alpha[K + 1];

	double f_filtr[K + 1];
	/*for (int k = 0; k <= K; k++) //Усредняем 
	{
		if (k > M && k < K - M)
		{
			f_filtr[k] += f_shum[k] * alpha[k];
		}
	}*/

	double w[10000];
	w[0] = 0;//Критерий зашумлённости, омега
	/*for (int k = 1; k <= N; k++)
	{
		w += abs(f_filtr[k] - f_filtr[k - 1]);
	}*/

	double d[10000];
	d[0] = 0;//Критерий близости, дельта
	/*for (int k = 0; k <= K; k++)
	{
		d += abs(f_filtr[k] - f_shum[k]);
	}
	d /= K;*/

	double lambda[L + 1]; //лямбды
	for (int l = 0; l <= 10; l++)
	{
		lambda[l] = l / L;
	}

	double J[11]; //Линейная свёртка критериев
	for (int i = 0; i < N; i++)
	{
		for (int k = K / 2; k <= K; k++) // функция генерации 2ой половины весов (50 - 100)
		{
			if (k != K / 2)
				alpha[k] = 0.5*random(0, 1 - alpha[k - 1]);
			else
				alpha[k] = random(0, 1);
		}

		for (int k = 0; k < K / 2; k++) // заполнение 1ой половины (0 - 100)
			alpha[k] = alpha[K - k];

		for (int k = 0; k <= K; k++) //Усредняем 
		{
			if (k > M && k < K - M)
			{
				f_filtr[k] += f_shum[k] * alpha[k];
			}
		}

		for (int k = 1; k <= N; k++)
		{
			w[i] += abs(f_filtr[k] - f_filtr[k - 1]);
		}

		for (int k = 0; k <= K; k++)
		{
			d[i] += abs(f_filtr[k] - f_shum[k]);
		}
		d[i] /= K;
	}

	auto _res_w = min_element(w, w+N);
	auto _res_d = min_element(d, d+N);

	double res_w = *_res_w;
	double res_d = *_res_d;

	for (int l = 0; l <= 10; l++)
	{
		J[l] = lambda[l] * res_w + (1 - lambda[l])*res_d;
	}

	double dist = abs(res_w) + abs(res_d);
}
