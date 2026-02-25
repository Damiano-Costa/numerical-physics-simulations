/*****************************************************************************/
/*
 * Questo programma simula un Gas Reticolare (Lattice Gas) su una griglia
 * bidimensionale di dimensione (GRID+1)*(GRID+1).
 * Le particelle si muovono secondo una logica di "Random Walk" con esclusione
 * di volume (non possono occupare la stessa cella).
 *
 * Il codice è strutturato in diverse sezioni controllate da direttive al
 * preprocessore (#ifdef):
 * 1. DTR: Calcola il coefficiente di diffusione D in funzione del tempo e
 * successivamente in funzione della densità (Rho). Salva i dati su file.
 * 2. EFFTAGLIA: Analizza come varia D al variare della dimensione della
 * griglia (GRID) per verificare gli effetti di taglia finita.
 * 3. TLIM: Verifica il comportamento asintotico per tempi lunghi.
 * 4. STATICFLUCT: Calcola le fluttuazioni statistiche ripetendo la
 * simulazione per molte "storie" diverse.
 * 
 * Il codice è pensato per eseguire una sola direttiva per volta
 * 
 * La funzione 'filler' posiziona casualmente le particelle all'inizio.
 * La funzione 'Event' esegue la dinamica (tentativi di spostamento).
 * La funzione 'DR' calcola lo spostamento quadratico medio per ottenere D.
 */
/*****************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stddef.h>

//Funzione che posiziona casualmente Npart particelle nel reticolo (senza sovrapposizioni)
//e inizializza i vettori delle posizioni
void filler(int **Lat, int **Peripos, int **TruePos, int GRID, int Npart, int **P0);

//Funzione che esegue un passo temporale della simulazione: prova a muovere le particelle
//in una direzione casuale e aggiorna le posizioni se la cella di arrivo è libera
void Event(int **Lat, int **Peripos, int **TruePos, int GRID, int Npart);

//Funzione che calcola lo spostamento quadratico medio
//rispetto alle posizioni iniziali
double DR(int **Pos0, int **Truepos, int Npart);

//Funzione per il calcolo della media dei coefficienti di diffusione su più storie
double mean(double *DRF, int Stories);

//Funzione per il calcolo della deviazione standard dei risultati
double SDT(double *DRF, int Stories, double M);


int main()
{
    srand(time(NULL));

    //Variabili generali
    int GRID = 100;                  //Dimensione del reticolo
    int maxpart = (GRID+1)*(GRID+1); //Numero massimo teorico di particelle (reticolo pieno)
    int Npart = 0.8*maxpart;         //Numero di particelle attuali (inizializzato a densità 0.8)
    
    //Allocazione memoria per gli array di posizione
    int **Pos0 = calloc(maxpart, sizeof(int*));     //Array per le posizioni al tempo t=0
    for(int i=0; i<maxpart; i++) {
        Pos0[i] = calloc(2, sizeof(int));
    }

    int **Lat = calloc(GRID+1, sizeof(int*));       //Matrice che rappresenta la griglia occupata
    for(int i=0; i<=GRID; i++) {
        Lat[i] = calloc(GRID+1, sizeof(int));
    }

    int **Truepos = calloc(maxpart, sizeof(int*));  //Array posizioni reali ( per calcolo D)
    for(int i=0; i<maxpart; i++) {
        Truepos[i] = calloc(2, sizeof(int));
    }

    int **Peripos = calloc(maxpart, sizeof(int*));  //Array posizioni periodiche (condizioni periodiche al bordo)
    for(int i=0; i<maxpart; i++) {
        Peripos[i] = calloc(2, sizeof(int));
    }

    double R2;                      //Variabile per lo spostamento quadratico medio
    double D;                       //Variabile per il coefficiente di diffusione
    FILE *file1, *file2;            //Puntatori ai file di output

/*
 * SEZIONE DTR: Calcolo Diffusione vs Tempo e Diffusione vs Densità
 */
#ifdef DTR
    //Parte 1: Diffusione in funzione del tempo
    file1 = fopen("CoeDiff.txt", "w");
    filler(Lat, Peripos, Truepos, GRID, Npart, Pos0); //Inizializzazione sistema
    
    for(int t=1; t<=10000; t++)     //Ciclo temporale
    {
        Event(Lat, Peripos, Truepos, GRID, Npart);
        R2 = DR(Pos0, Truepos, Npart);
        D = R2 / (4.0 * (double)t);
        fprintf(file1, "%i %f\n", t, D);
    }
    fclose(file1);

    //Parte 2: Diffusione in funzione della densità (Rho)
    file2 = fopen("CoediffvsRho.txt", "w");
    double Rhomin = 0.20;           //Densità minima
    double Rhomax = 0.90;           //Densità massima
    double step = 0.01;             //Passo di incremento densità
    int K = (int)((Rhomax - Rhomin) / step);
    double Rho = Rhomin;

    //Reset della griglia
    for(int k=0; k<=GRID; k++) {
        for(int i=0; i<=GRID; i++) {
            Lat[k][i] = 0;
        }
    }

    for(int j=0; j<=K; j++)         //Ciclo sulle diverse densità
    {
        Npart = Rho * (GRID+1) * (GRID+1);
        filler(Lat, Peripos, Truepos, GRID, Npart, Pos0);
        double time = 0;
        
        for(int t=1; t<=10000; t++) //Evoluzione temporale per ogni densità
        {
            Event(Lat, Peripos, Truepos, GRID, Npart);
            time = t;
        }
        
        R2 = DR(Pos0, Truepos, Npart);
        D = R2 / (4.0 * time);
        fprintf(file2, "%f %f\n", Rho, D);

        //Reset array per la prossima iterazione
        for(int k=0; k<=GRID; k++) {
            for(int i=0; i<=GRID; i++) {
                Lat[k][i] = 0;
            }
        }
        Rho += step;
    }
    fclose(file2);
#endif

/*
 * SEZIONE EFFTAGLIA: Analisi degli effetti di taglia finita del reticolo
 */
#ifdef EFFTAGLIA
    //Liberazione della memoria allocata inizialmente nel main
    for(int i=0; i<=GRID; i++)
    {
        free(Lat[i]);
    }
    free(Lat);
    for(int i=0; i<maxpart; i++)
    {
        free(Pos0[i]);
        free(Truepos[i]);
        free(Peripos[i]);
    }
    free(Pos0);
    free(Truepos);
    free(Peripos);

    int Stor = 5;                  //Numero di storie su cui mediare
    double *Dfin = calloc(Stor, sizeof(double));
    FILE *file5 = fopen("efftaglia.txt", "w");
    
    double Rho_eff = 0.5;           //Densità fissata
    double Tmax = 2000.0;           //Tempo di simulazione per ogni storia
    int L0 = 10;                    //Dimensione iniziale
    int Lmax = 50;                  //Numero di incrementi di taglia
    int *L = calloc(Lmax, sizeof(int)); //Array delle dimensioni L
    
    for(int i=0; i<Lmax; i++)
    {
        L[i] = L0;
        L0 += 10;
    }

    //Ciclo sulle diverse dimensioni L della griglia
    for(int k=0; k<Lmax; k++)
    {
        GRID = L[k];
        Npart = Rho_eff * (GRID+1) * (GRID+1);
        maxpart = Npart + 1;
        
        //Allocazione memoria per la nuova dimensione GRID
        Pos0 = calloc(maxpart, sizeof(int*));
        for(int i=0; i<maxpart; i++)
        {
            Pos0[i] = calloc(2, sizeof(int));
        }
        
        Lat = calloc(GRID+1, sizeof(int*));
        for(int i=0; i<=GRID; i++)
        {
            Lat[i] = calloc(GRID+1, sizeof(int));
        }
        
        Truepos = calloc(maxpart, sizeof(int*));
        for(int i=0; i<maxpart; i++)
        {
            Truepos[i] = calloc(2, sizeof(int));
        }
        
        Peripos = calloc(maxpart, sizeof(int*));
        for(int i=0; i<maxpart; i++)
        {
            Peripos[i] = calloc(2, sizeof(int));
        }

        //Ciclo sulle storie per calcolare la media
        for(int j=0; j<Stor; j++)
        {
            filler(Lat, Peripos, Truepos, GRID, Npart, Pos0);
            
            //Evoluzione temporale
            for(int t=1; t<=(int)Tmax; t++)
            {
                Event(Lat, Peripos, Truepos, GRID, Npart);
            }
            
            R2 = DR(Pos0, Truepos, Npart);
            D = R2 / (4.0 * Tmax);
            Dfin[j] = D;

            //Reset della matrice Lat per la prossima storia
            for(int q=0; q<=GRID; q++)
            {
                for(int z=0; z<=GRID; z++)
                {
                    Lat[q][z] = 0;
                }
            }
        }
        
        D = mean(Dfin, Stor);
        fprintf(file5, "%i %f\n", L[k], D);

        //Liberazione memoria alla fine del ciclo sulla taglia
        for(int i=0; i<=GRID; i++)
        {
            free(Lat[i]);
        }
        free(Lat);
        for(int i=0; i<maxpart; i++)
        {
            free(Pos0[i]);
            free(Truepos[i]);
            free(Peripos[i]);
        }
        free(Pos0);
        free(Truepos);
        free(Peripos);
    }
    
    free(L);
    free(Dfin);
    fclose(file5);
#endif

/*
 * SEZIONE TLIM: Verifica del limite asintotico temporale
 */
#ifdef TLIM
    FILE *file3 = fopen("Tlim.txt", "w");
    double tmax_lim = 1000;
    
    for(int j=1; j<4; j++)          //Ciclo per variare l'ordine di grandezza di tmax
    {
        filler(Lat, Peripos, Truepos, GRID, Npart, Pos0);
        
        for(int t=1; t<=tmax_lim; t++)
        {
            Event(Lat, Peripos, Truepos, GRID, Npart);
            R2 = DR(Pos0, Truepos, Npart);
            D = R2 / (4.0 * (double)t);
            fprintf(file3, "%i %f %i\n", t, D, j);
        }
        
        tmax_lim = 10 * tmax_lim;   //Incremento ordine di grandezza
        
        //Reset Griglia
        for(int k=0; k<=GRID; k++) {
            for(int i=0; i<=GRID; i++) Lat[k][i] = 0;
        }
    }
#endif

/*
 * SEZIONE STATICFLUCT: Analisi delle fluttuazioni statistiche
 */
#ifdef STATICFLUCT
    FILE *file4 = fopen("staticfluct.txt", "w");
    
    //Deallocazione e riallocazione per miniaturizzare il sistema (GRID=30)
    for(int i=0; i<=GRID; i++) free(Lat[i]);
    free(Lat);
    for(int i=0; i<maxpart; i++) {
        free(Pos0[i]); free(Truepos[i]); free(Peripos[i]);
    }
    free(Pos0); free(Truepos); free(Peripos);

    GRID = 30;                      //Nuova dimensione ridotta
    Npart = 0.5 * (GRID+1) * (GRID+1);
    maxpart = Npart + 1;

    //Riallocazione variabili ridotte
    Pos0 = calloc(maxpart, sizeof(int*));
    for(int i=0; i<maxpart; i++) Pos0[i] = calloc(2, sizeof(int));
    
    Lat = calloc(GRID+1, sizeof(int*));
    for(int i=0; i<=GRID; i++) Lat[i] = calloc(GRID+1, sizeof(int));
    
    Truepos = calloc(maxpart, sizeof(int*));
    for(int i=0; i<maxpart; i++) Truepos[i] = calloc(2, sizeof(int));
    
    Peripos = calloc(maxpart, sizeof(int*));
    for(int i=0; i<maxpart; i++) Peripos[i] = calloc(2, sizeof(int));

    int Stories = 1000;             //Numero di realizzazioni
    int Tmax_stat = 10000;
    double M_stat;                  //Media
    double sig;                     //Deviazione Standard
    double *DRF_stat = calloc(Stories, sizeof(double));

    for(int i=0; i<Stories; i++)
    {
        filler(Lat, Peripos, Truepos, GRID, Npart, Pos0);
        for(int t=1; t<=10000; t++)
        {
            Event(Lat, Peripos, Truepos, GRID, Npart);
        }
        R2 = DR(Pos0, Truepos, Npart);
        D = R2 / (4.0 * (double)Tmax_stat);
        DRF_stat[i] = D;
        fprintf(file4, "%i %f\n", i, D);

        //Reset Griglia
        for(int k=0; k<=GRID; k++) {
            for(int w=0; w<=GRID; w++) Lat[k][w] = 0;
        }
    }
    
    M_stat = mean(DRF_stat, Stories);
    sig = SDT(DRF_stat, Stories, M_stat);
    
    fclose(file4);
#endif
    
    return 0;
}

// Funzione che calcola la media aritmetica dei valori contenuti nell'array DRF
double mean(double *DRF, int Stories)
{
    double M = 0;               // Accumulatore per la somma
    int i;                      // Indice di iterazione

    for(i = 0; i < Stories; i++)
    {
        M += DRF[i];
    }

    return M / Stories;         // Restituisce il valore medio
}

// Funzione che calcola la deviazione standard campionaria dei valori
double SDT(double *DRF, int Stories, double M)
{
    double V = 0;               // Accumulatore per la varianza
    int i;                      // Indice di iterazione

    for(i = 0; i < Stories; i++)
    {
        // Somma dei quadrati degli scarti dalla media
        V += (M - DRF[i]) * (M - DRF[i]);
    }

    // Calcolo deviazione standard (normalizzazione su N-1 per campione statistico)
    V = sqrt(V / (Stories - 1));

    return V;
}

// Funzione che calcola lo Spostamento Quadratico Medio (MSD) di tutte le particelle
// utilizzando le coordinate reali rispetto a quelle iniziali
double DR(int **Pos0, int **Truepos, int Npart)
{
    double R = 0;               // Accumulatore per la somma degli spostamenti quadratici
    int i;                      // Indice di iterazione

    for(i = 0; i < Npart; i++)
    {
        // Calcolo della distanza quadrata R^2 = dx^2 + dy^2 per la particella i-esima
        R += (Truepos[i][0] - Pos0[i][0]) * (Truepos[i][0] - Pos0[i][0]) +
             (Truepos[i][1] - Pos0[i][1]) * (Truepos[i][1] - Pos0[i][1]);
    }

    return R / ((double)Npart); // Restituisce la media degli spostamenti quadratici
}

// Funzione che posiziona casualmente le particelle sulla griglia evitando sovrapposizioni.
// Inizializza tutte le strutture dati per le coordinate (Lat, Peripos, TruePos, P0)
void filler(int **Lat, int **Peripos, int **TruePos, int GRID, int Npart, int **P0)
{
    int j;                      // Indice della particella
    int x, y;                   // Coordinate temporanee

    for(j = 1; j <= Npart; j++)
    {
        // Generazione coordinate casuali con controllo di esclusione (cella vuota)
        do
        {
            x = rand() % (GRID + 1);
            y = rand() % (GRID + 1);
        } while(Lat[x][y] != 0); // Ripete se la cella è già occupata

        // Assegnazione della particella alla griglia e salvataggio coordinate
        Lat[x][y] = j;

        // Inizializzazione vettori posizioni
        Peripos[j-1][0] = x;    // Posizione periodica (nel box)
        Peripos[j-1][1] = y;
        TruePos[j-1][0] = x;    // Posizione reale
        TruePos[j-1][1] = y;
        P0[j-1][0] = x;         // Posizione iniziale (t=0)
        P0[j-1][1] = y;
    }
}

// Funzione che esegue un Passo: tenta di muovere ogni particella
// in una direzione casuale, rispettando le condizioni al contorno periodiche
// e il principio di esclusione di volume
void Event(int **Lat, int **PeriPos, int **Truepos, int GRID, int Npart)
{
    int DIM = GRID + 1;         // Dimensione effettiva dell'array
    int p;                      // Indice particella
    int i, j;                   // Coordinate attuali
    int n_i, n_j;               // Coordinate del vicino (proposta)
    int di, dj;                 // Spostamento lungo x e y
    int q;                      // Variabile per la scelta casuale della direzione

    for(p = 1; p <= Npart; p++)
    {
        // Recupero posizione attuale della particella p-esima
        i = PeriPos[p-1][0];
        j = PeriPos[p-1][1];

        // Scelta casuale della direzione di movimento
        q = (rand() % 4) + 1;
        di = 0;
        dj = 0;

        if (q == 1) di = 1;      // destra
        if (q == 2) dj = 1;      // alto
        if (q == 3) di = -1;     // sinistra
        if (q == 4) dj = -1;     // giù

        // Calcolo nuove coordinate con condizioni al contorno periodiche 
        n_i = (i + di + DIM) % DIM;
        n_j = (j + dj + DIM) % DIM;

        // Controllo se il sito di destinazione è libero
        if (Lat[n_i][n_j] == 0)
        {
            // Aggiornamento Griglia: libera vecchia cella, occupa nuova
            Lat[i][j] = 0;
            Lat[n_i][n_j] = p;

            // Aggiornamento coordinate periodiche
            PeriPos[p-1][0] = n_i;
            PeriPos[p-1][1] = n_j;

            // Aggiornamento coordinate reali
            Truepos[p-1][0] += di;
            Truepos[p-1][1] += dj;
        }
    }
}
