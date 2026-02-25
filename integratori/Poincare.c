/*****************************************************************************/
/*
 * Questo programma integra numericamente l'equazione del moto di un 
 * oscillatore forzato e smorzato tramite il metodo Runge-Kutta al quarto 
 * ordine (RK4).
 * Viene preso in input un indice per scegliere l'ampiezza della forzante.
 * Il programma genera due file: 'exem1.txt' per la sezione di Poincaré 
 * e 'exem2.txt' per la traiettoria continua nello 
 * spazio delle fasi.
 * La struct SF contiene lo stato del sistema e i parametri fisici.
 */
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//Struct per le variabili di stato e i parametri
struct SF
{
    double x;       //Posizione
    double v;       //Velocità
    double Om2;     //Omega al quadrato
    double Ttot;    //Tempo totale
    double F0;      //Ampiezza forzante
    double Gm;      //Smorzamento
    double Omm;     //omega forzante
};

//Integrazione numerica RK4
void RK4(struct SF *p, double Dt, double t);

//Calcolo dell'accelerazione
double Fi(struct SF *p, double t, double xn, double vn);

int main(int argc, char **argv)
{
    //Controllo argomenti
    if (argc != 2)
    {
        printf("errore negli argomenti\n");
        return 0;
    }

    if (argc == 2)
    {
        //Controllo validità indice input
        if (atoi(argv[1]) < 0 || atoi(argv[1]) > 3)
        {
            printf("errore nel valore degli argomenti");
            return 1;
        }

        if (atoi(argv[1]) >= 0 && atoi(argv[1]) <= 3)
        {
            //Variabili
            double t;                   //Tempo
            double Dt;                  //Passo temporale
            double r;                   
            FILE *output1, *output2;    //File di output

            output1 = fopen("exem1.txt", "w");
            output2 = fopen("exem2.txt", "w");

            struct SF p = {M_PI/2, 0, 1, 800, 0.90, 0.5, 2.0/3.0};
            
            //Parametri temporali legati alla forzante
            double TF0 = 2*M_PI/p.Omm;
            double Forzanti[4] = {1.07, 1.15, 1.47, 1.50};
            
            Dt = TF0/1000;
            int NDt = (int)(p.Ttot/Dt);
            double tol = Dt/4;
            
            t = 0;
            p.F0 = Forzanti[atoi(argv[1])];
            
            //Salvataggio condizioni iniziali
            double x0 = p.x;
            double v0 = p.v;            
            for (int j = 0; j < 10*NDt; j++)
            {
                t = j*Dt;
                RK4(&p, Dt, t);

                //Si ignorano i primi 20 periodi (transiente)
                if (t >= 20*TF0)
                {
                    //limito posizione in [-PI, PI]
                    while (p.x <= -M_PI) p.x = p.x + 2*M_PI;
                    while (p.x >= M_PI) p.x = p.x - 2*M_PI;
                    r = fmod(t, TF0);
                    //Salvataggio punti multipli del periodo
                    if (r < tol || TF0 - r < tol)
                    {
                        fprintf(output1, "%f %f\n", p.x, p.v);
                    }
                }
            }
            //Reset condizioni iniziali
            t = 0;
            p.x = x0;
            p.v = v0;
            //Ciclo per la traiettoria completa
            for (int i = 0; i < NDt; i++)
            {
                t = i*Dt;
                RK4(&p, Dt, t);
                //Si ignorano i primi 20 periodi
                if (t >= 20*TF0)
                {
                    while (p.x <= -M_PI) p.x = p.x + 2*M_PI;
                    while (p.x >= M_PI) p.x = p.x - 2*M_PI;
                    fprintf(output2, "%f %f\n", p.x, p.v);
                }
            }
            fclose(output1);
            fclose(output2);
        }
    }
}
//Funzione che calcola l'accelerazione del sistema
double Fi(struct SF *p, double t, double xn, double vn)
{
    return -p->Om2*sin(xn) - p->Gm*vn + p->F0*cos(p->Omm*t);
}

// Runge-Kutta
void RK4(struct SF *p, double Dt, double t)
{
    double xn, vn, Dx1, Dv1, Dx2, Dv2, Dx3, Dv3, Dx4, Dv4, ts;
    ts = t;
    xn = p->x;
    vn = p->v;
    //Calcolo incrementi
    Dx1 = vn*Dt;
    Dv1 = Fi(p, ts, xn, vn)*Dt;

    Dx2 = vn*Dt + Dv1*Dt*0.5;
    Dv2 = Fi(p, ts + 0.5*Dt, xn + 0.5*Dx1, vn + 0.5*Dv1)*Dt;

    Dx3 = vn*Dt + Dv2*Dt*0.5;
    Dv3 = Fi(p, ts + 0.5*Dt, xn + 0.5*Dx2, vn + 0.5*Dv2)*Dt;

    Dx4 = vn*Dt + Dv3*Dt;
    Dv4 = Fi(p, ts + Dt, xn + Dx3, vn + Dv3)*Dt;

    //Aggiornamento stato
    p->x = xn + (1.0/6.0)*Dx1 + (1.0/3.0)*Dx2 + (1.0/3.0)*Dx3 + (1.0/6.0)*Dx4;
    p->v = vn + (1.0/6.0)*Dv1 + (1.0/3.0)*Dv2 + (1.0/3.0)*Dv3 + (1.0/6.0)*Dv4;
}
