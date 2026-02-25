/*****************************************************************************/
/*
 * Questo codice simula un Random Walk bidimensionale (RW 2D) su reticolo.
 * L'output dipende dalle direttive del preprocessore:
 * - GRAPHIC: Calcola una singola traiettoria (x,y) per visualizzazione.
 * - FIT: Calcola lo spostamento quadratico medio <R^2> in funzione del tempo.
 * - HISTO: Costruisce un istogramma 2D delle posizioni a un tempo fissato.
 * Output su file "RW2d.txt".
 */
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

struct PQ 
{
    int x;          // Coordinata x
    int y;          // Coordinata y
    double d;       // Distanza dall'origine
};

void SM(struct PQ *s);
void Dist(struct PQ *s);

int main()
{
    struct PQ s = {0, 0, 0};
    srand(time(NULL));
    FILE *file;
    int Nstep = 300;
    int Walkers = 100000;
    
    file = fopen("RW2d.txt", "w");
    
    /************************** Modalità GRAPHIC **************************/
    /* Visualizzazione di una singola traiettoria traslata */
#ifdef GRAPHIC
    int gl = (int)((Nstep + 1) / 2); // Shift per centrare il grafico
    int x = s.x + gl;
    int y = s.y + gl;
    fprintf(file, "%i %i\n", x, y);
    
    for(int k = 0; k < Nstep + 1; k++)
    {
        SM(&s);
        x = s.x + gl;
        y = s.y + gl;
        fprintf(file, "%i %i\n", x, y);
    }
#endif

    /**************************** Modalità FIT ****************************/
    /* Calcolo di <R^2> vs t mediando su molti camminatori */
#ifdef FIT
    for(int t = 0; t < Nstep; t++)
    {
        double av = 0;
        // media  per ogni istante t
        for(int k = 0; k < Walkers; k++)
        {
            s.x = 0;
            s.y = 0;
            // Evoluzione fino al tempo t corrente
            for(int j = 0; j < t; j++)
            {
                SM(&s);
            }
            Dist(&s);
            av += s.d * s.d; // Somma quadratica
        }
        av = av / Walkers;
        fprintf(file, "%i %f\n", t, av);
    }
#endif

    /*************************** Modalità HISTO ***************************/
    /* Distribuzione spaziale 3D a tempo fissato */
#ifdef HISTO
    int TIMEHISTO = 120;
    int DimBin = (int)(3 * sqrt(TIMEHISTO)); // Range prendendo 3 sigma
    int Q = (2 * DimBin + 1);
    
    // Allocazione dinamica matrice 2D per i bin
    int **Bin = calloc(Q, sizeof(int*));
    for(int i = 0; i < Q; i++) 
    {
        Bin[i] = calloc(Q, sizeof(int));
    }
    
    // Simulazione 
    for(int k = 0; k < Walkers; k++)
    {
        s.x = 0;
        s.y = 0;
        for(int p = 0; p < TIMEHISTO; p++)
        {
            SM(&s); 
        }
        // Riempimento istogramma se dentro il range
        if(s.x >= -DimBin && s.x <= DimBin && s.y >= -DimBin && s.y <= DimBin)
        {
            Bin[s.x + DimBin][s.y + DimBin] += 1;
        }
    }
    
    // Scrittura matrice (formato compatibile con gnuplot splot)
    for(int i = 0; i < Q; i++)
    {
        for(int j = 0; j < Q; j++)
        {
            fprintf(file, "%i %i %i\n", i - DimBin, j - DimBin, Bin[i][j]);
        }
        fprintf(file, "\n"); // Riga vuota per separazione blocchi
    }
#endif

    fclose(file);
}

// Funzione Step Move: esegue un passo in una delle 4 direzioni
void SM(struct PQ *s)
{
    int q = rand() % 4;
    if (q == 0) s->x += 1;
    if (q == 1) s->x -= 1;
    if (q == 2) s->y += 1;
    if (q == 3) s->y -= 1;
}

// Calcola la distanza euclidea dall'origine
void Dist(struct PQ *s)
{
    s->d = sqrt(s->x * s->x + s->y * s->y); 
}
