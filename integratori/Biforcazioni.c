/*****************************************************************************/
/*
 * Diagramma di Biforcazione del Pendolo Smorzato e Forzato.
 *
 * Il programma analizza l'evoluzione del sistema al variare dell'ampiezza
 * della forzante (F0).
 * Input: Non riechiesto
 * Output: File "Biforcazioni.txt" con 4 colonne:
 * 1: Posizione (riportata tra -PI e PI)
 * 2: Velocità
 * 3: Ampiezza Forzante (Parametro di controllo)
 * 4: Indice k (per colorazione/distinzione traiettorie)
 */
/*****************************************************************************/

#include <math.h>
#include <stdio.h>

// Struttura parametri e stato
struct SF
{
	double x;       // Posizione
	double v;       // Velocità
	double Om2;     // Omega quadro 
	double Ttot;    // Tempo totale simulazione
	double F0;      // Ampiezza forzante (variabile nel loop)
	double Gm;      // Smorzamento
	double Omm;     // Frequenza forzante
};

// Funzione
void RK4(struct SF *p, double Dt, double t);
double Fi(struct SF *p, double t, double xn, double vn);

int main()
{
	// Inizializzazione parametri
	struct SF p = {M_PI / 2.0, 0, 1, 800, 0.90, 0.5, 2.0 / 3.0};
	// Calcolo passo temporale sincronizzato col periodo della forzante
	double TF0 = (2 * M_PI) / p.Omm;  // Periodo della forzante
	double Dt = TF0 / 100.0;          // 100 passi per periodo
	int NDt = (int)(p.Ttot / Dt);

	double t;
	
	// Parametri di variazione Forzante
	double F0start = 0.90;
	double F0end = 1.50;
	double F0step = 0.0001;

	FILE *output1 = fopen("Biforcazioni.txt", "w");

	// Loop esterno: varia le condizioni iniziali di velocità (k)
	for(int k = 0; k <= 10; k++)
	{
		// Loop intermedio: varia il parametro di controllo F0 (Biforcazione)
		for(p.F0 = F0start; p.F0 <= F0end; p.F0 += F0step)
		{
			// Reset condizioni iniziali per ogni nuovo F0
			p.x = M_PI / 2.0;
			p.v = k * (M_PI / 10.0);
			t = 0;

			// Loop interno: Integrazione temporale
			for(int j = 0; j < NDt; j++)
			{
				RK4(&p, Dt, t);
				t = (j + 1) * Dt; // Aggiornamento tempo preciso
				// j >= 2000: Salta il transiente iniziale.
				// j % 100 == 0: Poiché Dt = T/100, stampa esattamente una volta per periodo.
				if(j >= 2000 && (j % 100 == 0))
				{
					// limitazione dell'angolo tra -PI e PI
					while(p.x > M_PI)  p.x -= 2.0 * M_PI;
					while(p.x < -M_PI) p.x += 2.0 * M_PI;

					fprintf(output1, "%f %f %f %i\n", p.x, p.v, p.F0, k);
				}
			}
		}
	}

	fclose(output1);
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

	// aggiornamento finale
	p->x = xn + (1.0 / 6.0) * Dx1 + (1.0 / 3.0) * Dx2 + (1.0 / 3.0) * Dx3 + (1.0 / 6.0) * Dx4;
	p->v = vn + (1.0 / 6.0) * Dv1 + (1.0 / 3.0) * Dv2 + (1.0 / 3.0) * Dv3 + (1.0 / 6.0) * Dv4;
}
