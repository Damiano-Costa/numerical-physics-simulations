/*****************************************************************************/
/*
 * Simulazione dell'Oscillatore Armonico Semplice
 * Il programma richiede parametri iniziali e l'opzione del metodo di integrazione.
 * Calcola l'evoluzione temporale e l'energia meccanica, permettendo confronti
 * di stabilità e convergenza tra i diversi metodi d'integrazione.
 * * Opzioni:
 * 1-5: Integrazione standard (Eulero, Cromer, Punto Centrale, Verlet Auto, Verlet).
 * 6:   Analisi convergenza (Fit lineare dell'errore al variare di Dt) Utilizzare il comando grep per estrapolare dal txt
 *il metodo di cui si vuole fare il fit.
 */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// Struct per variabili di stato e parametri del sistema
struct SF {
    double x;       // Posizione corrente
    double v;       // Velocità corrente
    double Om2;     // Pulsazione al quadrato (k/m)
    double Ttot;    // Tempo totale
    double x0;      // Posizione step precedente (per Position Verlet)
};

//funzioni di integrazione
void Eulero(struct SF *p, double Dt);
void EuleroCrom(struct SF *p, double Dt);
void PuntoCentrale(struct SF *p, double Dt);
void VerletAuto(struct SF *p, double Dt);
void Verlet(struct SF *p, double Dt);

int main(int argc, char *argv[])
{
    // Controllo input
    if (argc != 6) {
        printf("Errore argomenti: <x0> <v0> <Om2> <Ttot> <Opzione>\n");
        return -1;
    }

    // Variabili simulazione
    double X0, V0, E0, Ek, t, DE, Dt;
    int k, NDt;
    int option = atoi(argv[5]);

    // Inizializzazione parametri da input
    struct SF p = { atof(argv[1]), atof(argv[2]), atof(argv[3]), atof(argv[4]) };

    // condizioni iniziali
    p.x0 = p.x;
    X0 = p.x; V0 = p.v;
    Dt = 0.01;
    t = 0;
    DE = 0;
    E0 = 0.5 * V0 * V0 + 0.5 * p.Om2 * X0 * X0; // Energia iniziale
    NDt = (int)(p.Ttot / Dt);
    
    // Array passi temporali per analisi errore (Opzione 6)
    double DT[5] = {0.005, 0.002, 0.001, 0.0005, 0.0001};

    /************************** Esecuzione Metodi  **************************/

    // Opzione 1: Eulero
    if (option == 1) {
        for (k = 0; k < NDt; k++) {
            Eulero(&p, Dt);
            printf("%f %f %f %f\n", p.x, p.v, t, DE);
            t += Dt;
            Ek = p.x * p.x * 0.5 * p.Om2 + p.v * p.v * 0.5;
            DE = Ek - E0;
        }
    }

    // Opzione 2: Eulero-Cromer
    if (option == 2) {
        for (k = 0; k < NDt; k++) {
            EuleroCrom(&p, Dt);
            printf("%f %f %f %f\n", p.x, p.v, t, DE);
            t += Dt;
            Ek = p.x * p.x * 0.5 + p.v * p.v * 0.5;
            DE = Ek - E0;
        }
    }

    // Opzione 3: Punto Centrale
    if (option == 3) {
        for (k = 0; k < NDt; k++) {
            PuntoCentrale(&p, Dt);
            printf("%f %f %f %f\n", p.x, p.v, t, DE);
            t += Dt;
            Ek = p.x * p.x * 0.5 + p.v * p.v * 0.5;
            DE = Ek - E0;
        }
    }

    // Opzione 4: Velocity Verlet
    if (option == 4) {
        for (k = 0; k < NDt; k++) {
            VerletAuto(&p, Dt);
            printf("%f %f %f %f\n", p.x, p.v, t, DE);
            t += Dt;
            Ek = p.x * p.x * 0.5 + p.v * p.v * 0.5;
            DE = Ek - E0;
        }
    }

    // Opzione 5: Verlet
    if (option == 5) {
        //calcolo passo indietro per Position Verlet
        double x = p.x, v = p.v;
        p.x0 = x - v * Dt;

        for (k = 0; k < NDt; k++) {
            Verlet(&p, Dt);
            printf("%f %f %f %f\n", p.x, p.v, t, DE);
            t += Dt;
            Ek = p.x * p.x * 0.5 + p.v * p.v * 0.5;
            DE = Ek - E0;
        }
    }

    // Opzione 6: studio della stabilità (Tutti i metodi su diversi Dt)
    if (option == 6) {
        // Test Eulero
        for (int j = 0; j < 5; j++)
	  {
            Dt = DT[j]; NDt = (int)(p.Ttot / Dt);
            for (int i = 0; i < NDt; i++)
	      {
		Eulero(&p, Dt);
		t += Dt;
	      }
            Ek = p.Om2 * p.x * p.x * 0.5 + p.v * p.v * 0.5;
            DE = Ek - E0;
            printf("Eul %f %f\n", Dt, DE);
            p.x = X0; p.v = V0; t = 0; DE = E0; // Reset
        }
	//ripristino le condizioni iniziali
        p.x = X0;
	p.v = V0;
	t = 0;
	DE = E0;

        // Test Eulero-Cromer
        for (int j = 0; j < 5; j++) {
            Dt = DT[j]; NDt = (int)(p.Ttot / Dt);
            for (int i = 0; i < NDt; i++)
	      {
	      EuleroCrom(&p, Dt);
	      t += Dt;
	      }
            Ek = p.Om2 * p.x * p.x * 0.5 + p.v * p.v * 0.5;
            DE = Ek - E0;
            printf("Crom %f %f\n", Dt, DE);
            p.x = X0; p.v = V0; t = 0; DE = E0;
        }
        p.x = X0;
	p.v = V0;
	t = 0;
	DE = E0;

        // Test Punto Centrale
        for (int j = 0; j < 5; j++) {
            Dt = DT[j]; NDt = (int)(p.Ttot / Dt);
            for (int i = 0; i < NDt; i++) { PuntoCentrale(&p, Dt); t += Dt; }
            Ek = p.Om2 * p.x * p.x * 0.5 + p.v * p.v * 0.5;
            DE = Ek - E0;
            printf("Punto-Centrale %f %f\n", Dt, DE);
            p.x = X0; p.v = V0; t = 0; DE = E0;
        }
        p.x = X0; p.v = V0; t = 0; DE = E0;

        // Test Velocity Verlet
        for (int j = 0; j < 5; j++) {
            Dt = DT[j]; NDt = (int)(p.Ttot / Dt);
            for (int i = 0; i < NDt; i++) { VerletAuto(&p, Dt); t += Dt; }
            Ek = p.Om2 * p.x * p.x * 0.5 + p.v * p.v * 0.5;
            DE = Ek - E0;
            printf("Ver-Auto %f %f\n", Dt, DE);
            p.x = X0; p.v = V0; t = 0; DE = E0;
        }
        p.x = X0;
	p.v = V0;
	t = 0;
	DE = E0;

        // Test Position Verlet
        for (int j = 0; j < 5; j++) {
            Dt = DT[j]; NDt = (int)(p.Ttot / Dt);
            p.x0 = p.x - p.v * Dt; // Bootstrap
            for (int i = 0; i < NDt; i++) { Verlet(&p, Dt); t += Dt; }
            Ek = p.Om2 * p.x * p.x * 0.5 + p.v * p.v * 0.5;
            DE = Ek - E0;
            printf("Verlet %f %f\n", Dt, DE);
            p.x = X0;
	    p.v = V0;
	    t = 0;
	    DE = E0;
        }
    }
    return 0;
}

/************************** Implementazione Metodi  **************************/

// Metodo di Eulero:
void Eulero(struct SF *p, double Dt) {
    double x = p->x;
    double v = p->v;
    p->x = x + v * Dt;
    p->v = v - p->Om2 * x * Dt;
}

// Metodo Eulero-Cromer: 
void EuleroCrom(struct SF *p, double Dt) {
    p->v = p->v - p->Om2 * p->x * Dt;
    p->x = p->x + p->v * Dt;
}

// Metodo Punto Centrale: 
void PuntoCentrale(struct SF *p, double Dt) {
    double v = p->v;
    p->v = v - p->Om2 * p->x * Dt;
    p->x = p->x + Dt * p->v * 0.5 + v * Dt * 0.5;
}

// Velocity Verlet:
void VerletAuto(struct SF *p, double Dt) {
    double xn = p->x;
    p->x = p->x + p->v * Dt - p->Om2 * p->x * 0.5 * Dt * Dt;
    double Fi12 = -1 * p->Om2 * xn - p->Om2 * p->x;
    p->v = p->v + Fi12 * 0.5 * Dt;
}

// Position Verlet: 
void Verlet(struct SF *p, double Dt) {
    double x_old = p->x0;
    double x_new = 2 * p->x - x_old - p->Om2 * p->x * Dt * Dt;
    p->v = (x_new - x_old) / (2 * Dt);
    p->x0 = p->x;
    p->x = x_new;
}
