/*****************************************************************************/
/*
 * Analisi dei Bacini di Attrazione (Pendolo Smorzato con Forzante).
 * Il programma genera una griglia di condizioni iniziali (X0, V0) e integra
 * il moto per classificare lo stato finale (con la variable controller).
 *
 * Input da riga di comando: Indice dell'array Forzante (0-3).
 * L'ordine delle forzanti è dalla più piccola proposta sul file alla più grande
 * Output: 3 colonne (X0, V0, controller) utilizzare la var controller per settare la palette di tipo Heatmap.
 * La variabile controller campiona la velocità finale se è positiva vale 1 e se è negativa o 0 è uguale 0.
 * Parametri:
 * NP: Numero di punti per lato della griglia (NP x NP simulazioni).
 * Controller: 1 se v_finale > 0, 0 se v_finale < 0.
 */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define NP 800           // Risoluzione griglia (quanti punti nell'effettivo)
#define Xin -(M_PI)      // Range iniziale X
#define Vin -M_PI        // Range iniziale V
#define Toll 0.001

struct SF
{
	double x;
	double v;
	double Om2;
	double Ttot;
	double F0;
	double Gm;
	double Omm;
};

// funzioni
void RK4(struct SF *p, double Dt, double t);
double Fi(struct SF *p, double t, double xn, double vn);


int main(int argc, char *argv[])
{
	// Controllo argomenti
	if(argc != 2)
	{
		printf("argomento invalido\n");
		return -1;
	}

	if(argc == 2)
	{
		// Allocazione array condizioni iniziali
		double X0[NP];
		double V0[NP];
		double Forzanti[4] = {1.07, 1.15, 1.47, 1.50};
		
		double t = 0;
		double Dt = 0.01;
		double passo = 2 * M_PI / (NP - 1);
		double xin = Xin;
		double vin = Vin;
		int controller;

		// Setup parametri fisici base
		struct SF p;
		p.Ttot = 100;
		p.Gm = 0.5;
		p.Omm = 2.0 / 3.0;
		p.Om2 = 1;
		
		int Nstep = (int)(p.Ttot / Dt);

		// Inizializzazione assi della griglia
		for(int j = 0; j < NP; j++)
		{
			X0[j] = xin;
			xin += passo;
		}

		for(int k = 0; k < NP; k++)
		{
			V0[k] = vin;
			vin += passo;
		}

		// Controllo validità input 
		int opzione = atoi(argv[1]);

		if(opzione > 3 || opzione < 0)
		{
			printf("errore\n");
		}
		
		if(opzione <= 3 && opzione >= 0)
		{
			// Aggiustamento tempo simulazione per casi specifici
			if(opzione > 1) 
			{
				p.Ttot = 93;
			}
			
			// Ricalcolo Nstep nel caso Ttot sia cambiato
			Nstep = (int)(p.Ttot / Dt);
			p.F0 = Forzanti[opzione];

			// Ciclo sulla griglia (Posizione iniziale i, Velocità iniziale j)
			for(int i = 0; i < NP; i++)
			{
				for(int j = 0; j < NP; j++)
				{
					// Reset condizioni singola simulazione
					p.x = X0[i];
					p.v = V0[j];
					t = 0;

					// Integrazione temporale
					for(int k = 0; k < Nstep + 1; k++)
					{
						RK4(&p, Dt, t);

						// All'ultimo step, valuta il "controller" e stampa
						if(k == Nstep)
						{
							if (p.v > 0) controller = 1;
							if (p.v < 0) controller = 0;
							
							printf("%f %f %i \n", X0[i], V0[j], controller);
						}
						
						t += Dt;
					}
				}
			}
		}
	}
	return 0;
}

/************************** Implementazione Funzioni **************************/

// Funzione Accelerazione
double Fi(struct SF *p, double t, double xn, double vn)
{
	double a;
	a = -p->Om2 * sin(xn) - p->Gm * vn + p->F0 * cos(p->Omm * t);
	return a;
}

// Metodo Runge-Kutta 4
void RK4(struct SF *p, double Dt, double t)
{
	double xn = p->x;
	double vn = p->v;

	// k1
	double k1x = vn;
	double k1v = Fi(p, t, xn, vn);

	// k2
	double k2x = vn + 0.5 * Dt * k1v;
	double k2v = Fi(p, t + 0.5 * Dt, xn + 0.5 * Dt * k1x, vn + 0.5 * Dt * k1v);

	// k3 
	double k3x = vn + 0.5 * Dt * k2v;
	double k3v = Fi(p, t + 0.5 * Dt, xn + 0.5 * Dt * k2x, vn + 0.5 * Dt * k2v);

	// k4 
	double k4x = vn + Dt * k3v;
	double k4v = Fi(p, t + Dt, xn + Dt * k3x, vn + Dt * k3v);

	// Aggiornamento stato
	p->x = xn + (Dt / 6.0) * (k1x + 2.0 * k2x + 2.0 * k3x + k4x);
	p->v = vn + (Dt / 6.0) * (k1v + 2.0 * k2v + 2.0 * k3v + k4v);
}


