#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include "lorentziana.h"
#include "savitzkygolay.h"

using namespace std;


struct tripla
{
	/*
	Estructura de datos que serve para almacenar os índices l, m e r dun triplete
	*/
	
	int l, m, r;
};

struct parametrosLorentz
{
	/*
	Estructura de datos que serve para almacenar os parámetros da curva lorentziana
	*/

	double l, A, w;
};

parametrosLorentz calculaparam(double *W, double y1, double y2, double y3, tripla p)
{
	/*
	Función para calcular os parámetros da lorentziana dentro do bucle correspondente
	*/

	double _w, beta, gamma, alpha, _l, _A, w12, w13, w23, y12, y23, y13, w1, w2, w3;
	parametrosLorentz param;

	w1 = *(W + p.l);
	w2 = *(W + p.m);
	w3 = *(W + p.r);

	w12 = w1 - w2;
	w13 = w1 - w3;
	w23 = w2 - w3;

	y12 = y1 - y2;
	y23 = y2 - y3;
	y13 = y1 - y3;



	_w = (pow(w1, 2) * y1 * y23 + pow(w3, 2) * y3 * y12 + pow(w2, 2) * y2 * (-y13)) / 
		(2 * w12 * y1 * y2 - 2 * y3 * (y1 * w13 - y2 * w23));

	_l = pow(-y12 / (y1 * pow(w1 - _w, 2) - y2 * pow(w2 - _w,2)), -0.5);

	_A = - 4 * w12 * w13 * w23 * y1 * y2 * y3 * (w1 * y1 * y23 + w3 * y12 * y3 + w2 * y2 * (- y13)) * _l /
		(pow(w12, 4) * pow(y1 * y2, 2) - 2 * pow(w12, 2) * y1 * y2 * (pow(w13, 2) * y1 + pow(w23, 2) * y2) *
		y3 + pow(pow(w13, 2) * y1 - pow(w23, 2) * y2, 2) * pow(y3, 2));


	param.l = _l;
	param.A = _A;
	param.w = _w;

	return param;
}

class espectro
{
	private:
		int np;	//número de elementos do espectro
		double *S;	//espectro
		double *S_dd;	//segunda derivada do espectro
		double delt; //parámetro umbral
		double *W;	//dominio
		double lamb; //parámetro das loretnzianas necesario para establecer a rexión libre de sinal

		double Score(tripla p) //función que computa a puntuación
		{
			double suma_L = 0, suma_R = 0;

			for(int i = p.l; i <= p.r; i++)
			{
				if(i <= p.m)
				{
					suma_L += abs(*(S_dd + i));
				}
				if(i >= p.m)
				{
					suma_R += abs(*(S_dd + i));
				}
			}
			if(suma_L < suma_R)
			{
				return suma_L;
			}
			else
			{
				return suma_R;
			}
		}
	
	public:
		vector<tripla> L; //tripletes sen filtrar
		vector<tripla> LL; //tripletes filtrados
		vector<parametrosLorentz> J; //vector que garda os parámetros de cada lorentziana 
		vector<double> scores;

		espectro(double *spectrum, double *derivada, double *dominio, int n, double delta)
		{
			/*
			Constructor
			*/

			np = n;
			S = spectrum;
			S_dd = derivada;
			W = dominio;
			delt = delta;
		}

		~espectro()
		{
			/*
			Destructor
			*/
			//delete[] S;
			//delete[] S_dd;
			L.clear();
			LL.clear();
		}

		void Guarda(double *spectrum, double *derivada, double *dominio, int n)
		{
			/*
			Función que almacena o espectro e computa a segunda derivada.
			*/

			np = n;
			S = spectrum;
			S_dd = derivada;
			W = dominio;

		}

		void Lee()
		{
			/*
			Función que almacena nun ficheiro os valores do espectro así como as súa derivada segunda
			*/

			ofstream ficheiroSaida;
			ficheiroSaida.open ("espectro.txt");
			ficheiroSaida << "Os valores que toma o espectro son: " << endl;
			for(int i = 0; i < np; i++)
			{
				ficheiroSaida << *(S + i) << endl;
			}
			ficheiroSaida << "Os valores que toma a segunda derivada son: " << endl;
			for(int i = 0; i < np; i++)
			{
				ficheiroSaida << *(S_dd + i) << endl;
			}
			ficheiroSaida.close();
		}

		void PlotPuntuacion()
		{
			/*
			Función que almacena nun ficheiro os valores necesarios para representar as puntuacións con gnuplot
			*/

			int k;
			ofstream ficheiro;
			ficheiro.open("plotPuntuacion.dat");
			for(int i = 0; i < L.size(); i++)
			{
				k = L[i].m;
				ficheiro << *(W + k) << "	" << scores[i] << endl;
			}
			ficheiro.close();
		}

		void PlotEspectro()
		{
			/*
			Función que almacena nun ficheiro os valores necesarios para representar o espectro con gnuplot
			*/

			ofstream ficheiroPlotEspectro;
			ficheiroPlotEspectro.open("plotEspectro.dat");
			for(int i = 0; i < np; i++)
			{
				ficheiroPlotEspectro << *(W + i) << "	" << setprecision(15) << *(S + i) << endl;
			}
			ficheiroPlotEspectro.close();
		}

		void PlotDerivada()
		{
			/*
			Función que almacena nun ficheiro os valores necesarios para representar a derivada do espectro
			*/
			ofstream ficheiroPlotDerivada;
			ficheiroPlotDerivada.open("plotDerivada.dat");
			for(int i = 0; i < np - 2; i++)
			{
				ficheiroPlotDerivada << *(W + i) << "	" << setprecision(15) << *(S_dd + i) << endl;
			}
			ficheiroPlotDerivada.close();
		}

		void PlotPosicionPicos()
		{
			/*
			Función que almacena nun ficheiro os valores necesarios para representar a posición dos picos
			*/

			ofstream ficheiroPicos;
			ficheiroPicos.open("plotPicos.dat");
			for(int i = 0; i < np; i++)
			{
				for(int j = 0; j < LL.size(); j++)
				{
					if(i == LL[j].m)
					{
						ficheiroPicos << *(W + i) << "	" << setprecision(15) << 0.04 << endl;
					}
				}
			}
			ficheiroPicos.close();	
		}

		void PlotNovoEspectro()
		{
			/*
			Función que almacena nun ficheiro os valores necesarios para representar o espectro reconstruido
			*/

			double *Q;
			Q = new double[np];


			for(int i = 0; i < np; i++)
			{
				for(int j = 0; j < J.size(); j++)
				{
					*(Q + i) += J[j].A * J[j].l / (pow(J[j].l, 2) + pow((*(W + i) - J[j].w), 2));
				}
			}

			for(int i = 0; i < J.size(); i++)
			{
				cout << J[i].A << "	" << J[i].l << "	" << J[i].w << endl;
			}

			ofstream ficheiroEspectro;
			ficheiroEspectro.open("plotNovoEspectro.dat");
			for(int i = 0; i < np; i++)
			{
				ficheiroEspectro << *(W + i) << "	" << setprecision(15) << *(Q + i) << endl;
			}
			ficheiroEspectro.close();	

			delete[] Q;
		}

		void Tripletes()
		{
			/*
			Función que calcula os tripletes do espectro
			*/

			int l, m, r;
			tripla p;
			double score, max_score = 0.0, score_media = 0.0;


			for(int i = 1; i < np - 3; i++)
			{
				if((*(S_dd + i) < 0 && *(S_dd + i - 1) > *(S_dd + i)) && *(S_dd + i + 1) > *(S_dd + i)) //buscamos mínimos locais negativos
				{
					m = i;

					for(int j = m - 1; j >= 1; j--)
					{
						if(((*(S_dd + j - 1) >= 0 && *(S_dd + j) < 0)) || ((*(S_dd + j - 1) <= *(S_dd + j) && *(S_dd + j + 1) < *(S_dd + j))))//raíz máis próxima ou máximo local máis próximo
						{
							l = j;
							break;
						}
					}

					for(int j = m + 1; j < np - 3; j++)
					{
						if(((*(S_dd + j + 1) >= 0 && *(S_dd + j) < 0)) || ((*(S_dd + j - 1) < *(S_dd + j) && *(S_dd + j) >= *(S_dd + j + 1))))
						{
							r = j;
							break;
						}
					}

					p.l = l;
					p.m = m;
					p.r = r; //almaceno en p l, m e r
					score = Score(p);

					L.push_back(p); //introduzo a tripla no vector
					scores.push_back(score);
				}
			}

			for(int i = 0; i < scores.size(); i++)
			{
				if(scores[i] > max_score)
				{
					max_score = scores[i];
				}
				score_media += scores[i];
			}

			for(int i = 0; i < scores.size(); i++)
			{
				scores[i] /= max_score;
			}
			score_media /= max_score * scores.size();

			double sdev = 0.0;
			int cont = 0;

			//Determinación da zona sen sinal:
			double rms;
			vector<double> RMS;

			for(int i = 0; i < np; i++)
			{
				rms += pow(*(S + i), 2);
				if(i % (np / 20) == 0)
				{
					rms = sqrt(rms);
					RMS.push_back(rms);
					rms = 0.0;
				}
			}

			sort(RMS.begin(), RMS.end());

			double mediana = RMS[RMS.size() / 2];

			//O de abaixo é para calcular a desviación típica das puntuacións dos tripletes en zonas sen ruído.
			//Empregando a mediana da desviación típica do ruído funciona perfectamente e automatizo este paso.

			/*
			for(int i = 0; i < L.size(); i++)
			{
				if(*(W + L[i].m) < sdev_r && *(W + L[i].m) > sdev_l) //igual sobra lambda dentro desta clase
				{
					sdev += pow((scores[i] - score_media), 2);
					cont += 1;
				}
			}

			sdev = sqrt(sdev / cont);*/
			
			for(int i = 0; i < scores.size(); i++)
			{
				if(scores[i] >= score_media + delt * 2 * mediana)
				{
					LL.push_back(L[i]);
				}
			}
		}

		void AproximacionProporcional()
		{
			/*
			Función que calcula os parámetros do espectro
			*/

			double Yl, Ym, Yr, yl, ym, yr, suml, summ, sumr;
			parametrosLorentz param;
			int maxIt = 1000, contit = 0;
			bool it;

			for(int i = 0; i < LL.size(); i++)
			{
				/*
				Neste bucle facemos unha primeira aproximación dos parámetros para comezar a iterar
				*/

				param = calculaparam(W, *(S + LL[i].l), *(S + LL[i].m), *(S + LL[i].r), LL[i]);

				J.push_back(param);
			}

			for(int i = 0; i < maxIt; i ++)
			{
				// Aquí iteramos ata que o resultado sexa o desexado

				for(int j = 0; j < LL.size(); j++)
				{
					it = true;

					Yl = *(S + LL[j].l);
					Ym = *(S + LL[j].m);
					Yr = *(S + LL[j].r);

					for(int k = 0; k < J.size(); k++)
					{
						yl = J[k].A * J[k].l / (pow(J[k].l,2) + pow((LL[j].l - J[k].w), 2));
						ym = J[k].A * J[k].l / (pow(J[k].l,2) + pow((LL[j].m - J[k].w), 2));
						yr = J[k].A * J[k].l / (pow(J[k].l,2) + pow((LL[j].r - J[k].w), 2));

						if(!(yl <= *(S + LL[j].l) && ym <= *(S + LL[j].m) && yr <= *(S + LL[j].r)))
						{
							it = false;
						}
					}

					if(it)
					{
						contit += 1;
						suml = 0;
						summ = 0;
						sumr = 0;
						for(int k = 0; k < J.size(); k++)
						{
							//	Calculo as sumas
							suml += J[k].A * J[k].l / (pow(J[k].l, 2) + pow((*(W + LL[j].l) - J[k].w), 2));
							summ += J[k].A * J[k].l / (pow(J[k].l, 2) + pow((*(W + LL[j].m) - J[k].w), 2));
							sumr += J[k].A * J[k].l / (pow(J[k].l, 2) + pow((*(W + LL[j].r) - J[k].w), 2));
						}

						//	Calculo as novas alturas
						Yl *= *(S + LL[j].l) / suml;
						Ym *= *(S + LL[j].m) / summ;
						Yr *= *(S + LL[j].r) / sumr;

						param = calculaparam(W, Yl, Ym, Yr, LL[j]);

						J[j] = param;
					}
				}
			}

			cout << "iteracions: " << contit / 100 << endl;
		}
};

int main()
{
	double *s, *dd, *w = 0;
	int num_w = 30000, num_picos = 100;
	double lambda = 0.0005, A = 1, w_max = 1, rho = 700;
	double delta = 3;
	bool ruido = true;
	/*
	s = new double[num_w];
	dd = new double[num_w];

	//xeración do espectro a partir da función "lorentz" pertencente ao encabezado "lorentziana.h"
	w = new double[num_w];

	for(int i = 1; i < num_w; i++)
	{
		*(w + i) = *w + i * w_max/num_w;
	}

	lorentz(s, w, lambda, A, rho, num_picos, num_w, w_max, ruido);

	sgconvolution(s, dd, 11, num_w, SG11);*/


	ifstream file ("PeakPicking.csv");

	double x, y;
	vector<double> X, Y;
	while(file >> x >> y)
	{
		X.push_back(x);
		Y.push_back(y);
	}
	file.close();

	num_w = X.size();


	vector<double> YY;
	//Engadimos ruído
	double desviaciontipica = (0.2 / 5) / rho;


	ofstream numpicos;
	numpicos.open("numpicos.dat");
	for(delta = 0; delta <= 25; delta += 0.1)
	{
		for(double sigma = 0.1; sigma >= 0.000000001; sigma /= 1.1)
		{



	for(int i =0; i < Y.size(); i++)
	{
		YY.push_back(gaussrand(sigma) + Y[i]);
	}

	s = &YY[0];
	w = &X[0];

	dd = new double[num_w];
	sgconvolution(s, dd, 5, num_w, SG5);

	espectro spec(s, dd, w, num_w, delta);

	//spec.Lee();
	spec.Tripletes();

	numpicos << delta << "		" << sigma << "		" << spec.LL.size() << endl;


	spec.L.clear();
	spec.LL.clear();
	YY.clear();
	cout << delta << endl;
}}

numpicos.close();
/*
	spec.PlotEspectro();
	spec.PlotDerivada();
	spec.PlotPuntuacion();
	spec.PlotPosicionPicos();
	spec.AproximacionProporcional();
	spec.PlotNovoEspectro();
	cout << spec.LL.size() << "	"<< spec.L.size() << "	" << spec.J.size();*/

	delete[] dd;
	X.clear();
	Y.clear();

	return 0;
}
