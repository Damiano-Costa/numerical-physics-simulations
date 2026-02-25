/*****************************************************************************/
/*
 * Questo programma simula un random walk (cammino casuale) unidimensionale
 * Il comportamento del programma è determinato da direttive al preprocessore
 * passate in fase di compilazione (es. -DGRAPH):
 * - GRAPH: Simula un singolo camminatore e salva la traiettoria (posizione 
 * vs tempo) su file.
 * - FIT: Calcola lo spostamento quadratico medio <x^2> su vari camminatori
 * al variare del tempo, per verificare l'andamento come t.
 * - HISTO: Costruisce la distribuzione di probabilità (istogramma) delle 
 * posizioni dei camminatori ad un tempo fissato.
 * I risultati sono salvati nel file "exem2.txt".
 */
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

//Struct per la gestione del passo casuale
struct PQ 
{
    int p;      //Valore del passo (+1 o -1)
    double q;   //Variabile per il numero casuale
};

//Funzione che genera il passo casuale
void SM(struct PQ *S);

int main()
{
    //Variabili generali
    int NS = 100;           //Numero di passi temporali
    int Walkers = 1000;     //Numero di camminatori (per le statistiche)
    struct PQ S;            //Struct per il passo
    FILE *output1;          //File di output

    //Apertura del file
    output1 = fopen("exem2.txt", "w");
    
    //Inizializzazione
    S.p = 0;
    srand(time(NULL));
    
    //Allocazione array posizioni
    int *Position = calloc(NS, sizeof(int));
    int p0 = S.p;

    /************************** Modalità GRAPH  **************************/
    /* Visualizzazione della traiettoria di un singolo camminatore */
#ifdef GRAPH
    //Generazione del cammino
    for (int k = 0; k < NS; k++)
    {
        Position[k] = p0;
        SM(&S);
        p0 = p0 + S.p;
    }

    //Scrittura su file della coppia (posizione, tempo)
    for (int j = 0; j < NS; j++)
    {
        fprintf(output1, "%i %i\n", Position[j], j);
    }
    
    //Pulizia memoria e chiusura
    free(Position);
    fclose(output1);
#endif

    /************************** Modalità FIT  ****************************/
    /* Calcolo dello spostamento quadratico medio <x^2> in funzione del tempo */
#ifdef FIT
    //Variabili specifiche per il fit
    double *Avw = calloc(Walkers, sizeof(double));
    int *TempFit = calloc(NS, sizeof(int));
    double avg = 0;

    //Inizializzazione array dei tempi
    for (int h = 0; h < NS; h++)
    {
        TempFit[h] = h;
    }

    //Ciclo sui tempi (per ogni istante temporale calcoliamo la media)
    for (int j = 0; j < NS; j++)
    {
        //Ciclo sui camminatori 
        for (int i = 0; i < Walkers; i++)
        {
            //Simulazione del cammino fino al tempo TempFit[j]
            for (int k = 0; k < TempFit[j]; k++)
            {
                Position[k] = p0;
                SM(&S);
                p0 += S.p;
            }
            //Accumulo del quadrato della distanza
            avg += p0 * p0;
            
            //Reset della posizione per il prossimo camminatore
            p0 = 0;
        }
        
        //Calcolo della media e scrittura su file
        avg = avg / Walkers;
        fprintf(output1, "%f %i\n", avg, TempFit[j]);
        avg = 0.0;
    }
    
    //Pulizia memoria e chiusura
    free(Avw);
    free(TempFit);
    free(Position);
    fclose(output1);
#endif

    /************************** Modalità HISTO  **************************/
    /* Calcolo dell'istogramma delle posizioni a un tempo fissato */
#ifdef HISTO
    //Ridefinizione parametri per statistica più robusta
    NS = 1000;
    Walkers = 100000;
    int temps = 850; //Tempo fissato per l'osservazione
    
    //Riallocazione memoria per il numero maggiore di passi
    int *Nstep = realloc(Position, NS * sizeof(int));
    Position = Nstep;

    //Definizione dei bin per l'istogramma
    int DimBin = (int)(3 * sqrt(temps)); //Range basato su tre sigma=3*sqrt(t)
    int Q = 2 * (DimBin + 1);
    int *ArrBin = calloc(Q, sizeof(int));
    
    //Array di appoggio per i tempi
    int *Temp = calloc(NS, sizeof(int));
    for (int h = 0; h < NS; h++)
    {
        Temp[h] = h;
    }
    //Ciclo sui camminatori
    for (int k = 0; k < Walkers; k++)
    {
        //Evoluzione temporale del singolo camminatore
        for (int j = 0; j < NS; j++)
        {
            Position[j] = p0;
            SM(&S);
            p0 += S.p;

            //Se siamo al tempo fissato 'temps', aggiorniamo l'istogramma
            if (Temp[j] == temps)
            {
                //Controllo se la particella è dentro il range del binning(3 sigma)
                if (Position[j] < (DimBin + 1) && Position[j] > -(DimBin + 1))
                {
                    //Shift dell'indice per usare array C (da 0 a Q)
                    Position[j] += DimBin;
                    ArrBin[Position[j]] += 1;
                }
            }
        }
        
        //Reset array posizioni (non strettamente necessario ma pulito)
        for (int f = 0; f < NS; f++)
        {
            Position[f] = 0;
        }
        
        //Reset variabili camminatore
        S.p = 0;
        p0 = S.p;
    }
   
    //Scrittura dell'istogramma su file
    for (int k = 0; k < Q; k++)
    {
        fprintf(output1, "%i %i\n", k, ArrBin[k]);
    }
    
    //Pulizia memoria e chiusura
    free(ArrBin);
    free(Temp);
    free(Position);
    fclose(output1);
#endif
}

//Funzione che genera il passo del random walk (+1 o -1) con probabilità 0.5
void SM(struct PQ *S)
{
    S->q = (double)rand() / RAND_MAX;
    
    if (S->q >= 0.5) S->p = 1;
    
    if (S->q < 0.5) S->p = -1; 
}
