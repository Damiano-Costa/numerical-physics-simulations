/*****************************************************************************/
/*
 * Simulazione Pendolo Smorzato con Forzante.
 * Risolve l'equazione differenziale completa:
 * a(t) = -Om^2*sin(x) - Gamma*v + F0*cos(Om_force * t)
 *
 * Input: x0, v0, Om^2, Ttot, F0, Gamma, Omega_forzante, Opzione.
 * Output: Stampa x, v, t, Energia (o Delta Energia).
 *
 * Opzioni:
 * 1: Eulero
 * 2: Eulero-Cromer
 * 3: Runge-Kutta 2
 * 4: Runge-Kutta 4
 * 5: Analisi convergenza (Fit DE vs Dt) per RK2 e RK4.
 */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// Struttura parametri e stato
struct SF
{
	double x;       // Posizione
	double v;       // Velocità
	double Om2;     // Omega quadro 
	double Ttot;    // Tempo totale
	double F0;      // Ampiezza forzante
	double Gm;      // Coefficiente smorzamento (Gamma)
	double Omm;     // Omega della forzante
};

// funzioni
double Fi(struct SF *p, double t, double xn, double vn);
void Eulero(struct SF *p, double Dt, double t);
void EuleroCrom(struct SF *p, double Dt, double t);
void RK2(struct SF *p, double Dt, double t);
void RK4(struct SF *p, double Dt, double t);

int main(int argc, char *argv[])
{
	// Controllo argomenti
	if(argc != 9)
	{
		printf("errore negli argomenti\n");
		return 0;
	}

	if(argc == 9)
	{
		double E0, Ek, t, DE, Dt;
		int k, NDt;
		int option = atoi(argv[8]);

		// Inizializzo struct
		struct SF p = {
			atof(argv[1]),
			atof(argv[2]),
			atof(argv[3]),
			atof(argv[4]),
			atof(argv[5]),
			atof(argv[6]),
			atof(argv[7])
		};

		// Array passi temporali (Opzione 5)
		double DT[5] = {0.2, 0.1, 0.05, 0.025, 0.0125};
		
		Dt = 0.01;
		t = 0;
		double X0 = p.x;
		double V0 = p.v;
		
		// Energia Iniziale
		E0 = 0.5 * p.v * p.v + p.Om2 * (1.0 - cos(p.x));
		NDt = (int)(p.Ttot / Dt);

		/************************* Esecuzione Metodi *************************/

		// Opzione 1: Eulero
		if(option == 1)
		{
			for(k = 0; k < NDt; k++)
			{
				Eulero(&p, Dt, t);
				Ek = 0.5 * p.v * p.v + p.Om2 * (1.0 - cos(p.x));
				DE = Ek - E0;
				printf("%f %f %f %f\n", p.x, p.v, t, DE);
				t = t + Dt;
			}
		}

		// Opzione 2: Eulero-Cromer
		if(option == 2)
		{
			for(k = 0; k < NDt; k++)
			{
				EuleroCrom(&p, Dt, t);
				Ek = 0.5 * p.v * p.v + p.Om2 * (1.0 - cos(p.x));
				DE = Ek - E0;
				printf("%f %f %f %f\n", p.x, p.v, t, DE);
				t = t + Dt;
			}
		}

		// Opzione 3: Runge-Kutta 2
		if(option == 3)
		{
			for(k = 0; k < NDt; k++)
			{
				RK2(&p, Dt, t);
				Ek = 0.5 * p.v * p.v + p.Om2 * (1.0 - cos(p.x));
				DE = Ek - E0;
				printf("%f %f %f %f\n", p.x, p.v, t, Ek);
				t = t + Dt;
			}
		}

		// Opzione 4: Runge-Kutta 4
		if(option == 4)
		{
			for(k = 0; k < NDt; k++)
			{
				RK4(&p, Dt, t);
				Ek = 0.5 * p.v * p.v + p.Om2 * (1.0 - cos(p.x));
				DE = Ek - E0;
				printf("%f %f %f %f\n", p.x, p.v, t, Ek);
				t = t + Dt;
			}
		}

		// Opzione 5:(Fit DE vs Dt)
		if(option == 5)
		{
			// Loop RK2
			for(int k = 0; k < 5; k++)
			{
				Dt = DT[k];
				NDt = (int)(p.Ttot / Dt);
				for(int j = 0; j < NDt; j++)
				{
					RK2(&p, Dt, t);
					t += Dt;
				}
				Ek = 0.5 * p.v * p.v + p.Om2 * (1.0 - cos(p.x));
				DE = Ek - E0;
				printf("RK2 %.15e %.15e\n", Dt, DE);  //scelto il formato %.15e poiche RK è preciso
				
				// Reset
				p.x = X0;
				p.v = V0;
				t = 0;
			}

			// Reset globale intermedio
			p.x = X0;
			p.v = V0;
			t = 0;

			// Loop convergenza RK4
			for(int k = 0; k < 5; k++)
			{
				Dt = DT[k];
				NDt = (int)(p.Ttot / Dt);
				for(int j = 0; j < NDt; j++)
				{
					RK4(&p, Dt, t);
					t += Dt;
				}
				Ek = 0.5 * p.v * p.v + p.Om2 * (1.0 - cos(p.x));
				DE = Ek - E0;
				printf("RK4 %.15e %.15e\n", Dt, DE);
				
				// Reset ()
				p.x = X0;
				p.v = V0;
				t = 0;
			}
		}
	}
	return 0;
}

/************************** Implementazione Metodi  **************************/

// Funzione accelerazione a(x, v, t)
double Fi(struct SF *p, double t, double xn, double vn)
{
	double a;
	a = -p->Om2 * sin(xn) - p->Gm * vn + p->F0 * cos(p->Omm * t);
	return a;
}

void Eulero(struct SF *p, double Dt, double t)
{
	double x, v;
	v = p->v;
	x = p->x;
	p->x = x + v * Dt;
	p->v = v - p->Om2 * sin(x) * Dt - p->Gm * v * Dt + p->F0 * cos(p->Omm * t) * Dt;
}

void EuleroCrom(struct SF *p, double Dt, double t)
{
	p->v = p->v - p->Om2 * sin(p->x) * Dt - p->Gm * p->v * Dt + p->F0 * cos(p->Omm * t) * Dt;
	p->x = p->x + p->v * Dt;
}

void RK2(struct SF *p, double Dt, double t)
{
	double xn, vn, Dx, Dv;
	xn = p->x;
	vn = p->v;
	
	// Step 1
	Dx = vn * Dt;
	Dv = Fi(p, t, xn, vn) * Dt;
	
	// Aggiornamento
	p->x = p->x + Dt * vn + 0.5 * Dt * Dv;
	p->v = p->v + Fi(p, t + 0.5 * Dt, xn + 0.5 * Dx, vn + 0.5 * Dv) * Dt;
}

void RK4(struct SF *p, double Dt, double t)
{
	double xn, vn, Dx1, Dv1, Dx2, Dv2, Dx3, Dv3, Dx4, Dv4, ts;
	ts = t;
	xn = p->x;
	vn = p->v;

	// K1
	Dx1 = vn * Dt;
	Dv1 = Fi(p, ts, xn, vn) * Dt;

	// K2
	Dx2 = vn * Dt + Dv1 * Dt * 0.5;
	Dv2 = Fi(p, ts + 0.5 * Dt, xn + 0.5 * Dx1, vn + 0.5 * Dv1) * Dt;

	// K3
	Dx3 = vn * Dt + Dv2 * Dt * 0.5;
	Dv3 = Fi(p, ts + 0.5 * Dt, xn + 0.5 * Dx2, vn + 0.5 * Dv2) * Dt;

	// K4
	Dx4 = vn * Dt + Dv3 * Dt;
	Dv4 = Fi(p, ts + Dt, xn + Dx3, vn + Dv3) * Dt;

	// aggiornamento
	p->x = xn + (1.0 / 6.0) * Dx1 + (1.0 / 3.0) * Dx2 + (1.0 / 3.0) * Dx3 + (1.0 / 6.0) * Dx4;
	p->v = vn + (1.0 / 6.0) * Dv1 + (1.0 / 3.0) * Dv2 + (1.0 / 3.0) * Dv3 + (1.0 / 6.0) * Dv4;
}
