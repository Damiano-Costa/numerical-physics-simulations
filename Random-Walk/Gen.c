/*****************************************************************************/
/*
 * Questo programma simula un Random Walk 1D confrontando 4 generatori 
 * pseudo-casuali (LEcuyer1/2, MinStand, InfRand).
 * Genera 3 file: graph.txt (traiettoria), fit.txt (MSD vs t) e histo.txt 
 * (distribuzione finale). Usare grep per estrarre la materia d'interesse
 * Richiede un seme intero da riga di comando.
 */
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

struct PQ {
    int p;          // Passo (+1 o -1)
    double q;       // Variabile casuale uniformemente distribuita
};

// Prototipi dei generatori
int LEcuyer1(int In);
int LEcuyer2(int In);
int MINSTAND(int In);
int INFRAND(int In);
void SM(struct PQ *S, int *seed, int method);

int main(int argc, char **argv)
{
    if(argc != 2) {
        printf("errore negli argomenti\n");
        return 0;
    }
    
    // Inizializzazione seme da input o default
    int I0 = atoi(argv[1]);
    if (I0 == 0) I0 = 123456789;

    struct PQ S;
    S.p = 0;
    FILE *filegraph = fopen("graph.txt", "w");
    FILE *filefit = fopen("fit.txt", "w");
    FILE *filehisto = fopen("histo.txt", "w");

    if (filegraph == NULL || filefit == NULL || filehisto == NULL) {
        printf("Errore apertura file\n");
        return 1;
    }

    char *NomiMetodi[4] = {"LEcuyer1", "LEcuyer2", "MINSTAND", "INFRAND"};
    
    // Ciclo principale: esegue la simulazione per tutti i 4 metodi
    for(int metodo = 0; metodo < 4; metodo++)
    {
        int seed = I0; // Reset del seed per confronto equo
        
        /************************** Traiettoria singola **************************/
        int NS = 1000;
        int p0 = 0;
        int *Position = calloc(NS, sizeof(int));
        
        // Calcolo evoluzione temporale singolo camminatore
        for (int k = 0; k < NS; k++) {
            Position[k] = p0;
            SM(&S, &seed, metodo);
            p0 = p0 + S.p;
        }
        // Scrittura su file (etichetta metodo, posizione, tempo)
        for(int j = 0; j < NS; j++) {
            fprintf(filegraph, "%s %i %i\n", NomiMetodi[metodo], Position[j], j);
        }
        free(Position);

        /******************************* Fit *********************************/
        /* Calcolo dello spostamento quadratico medio <x^2> */
        seed = I0;
        NS = 100;
        int Walkers = 1000;
        double *MSD = calloc(NS, sizeof(double));

        for(int i = 0; i < Walkers; i++) {
            int p = 0;
            for(int k = 0; k < NS; k++) {
                SM(&S, &seed, metodo);
                p += S.p;
                MSD[k] += (double)(p*p); // Accumulo x^2
            }
        }
        for(int j = 0; j < NS; j++) {
            fprintf(filefit, "%s %f %i\n", NomiMetodi[metodo], MSD[j]/Walkers, j);
        }
        free(MSD);
        
        /****************************** Istogramma *******************************/
        /* Distribuzione di probabilità al tempo fissato 'temps' */
        seed = I0;
        NS = 1000;
        Walkers = 100000; // Alta statistica per l'istogramma
        int temps = 850;
        
        // Definizione range istogramma basato su sigma ~ sqrt(t)
        int DimBin = (int)(3*sqrt(temps)); //prendo tutto ciò che cade entro 3 sigma
        int Q = 2*(DimBin+1);
        int *ArrBin = calloc(Q, sizeof(int));
        int *Temp = calloc(NS, sizeof(int));
        Position = calloc(NS, sizeof(int)); 

        for (int h = 0; h < NS; h++) Temp[h] = h;

        for(int k = 0; k < Walkers; k++) {
            p0 = 0;
            for(int j = 0; j < NS; j++) {
                SM(&S, &seed, metodo);
                p0 += S.p;
                Position[j] = p0;
                
                // Se siamo al tempo di osservazione
                if(Temp[j] == temps) {
                    // Check se la particella è nel range
                    if(Position[j] < (DimBin+1) && Position[j] > -(DimBin+1)) {
                        // Shift indice per array C (da 0 a Q)
                        int index = Position[j] + DimBin;
                        if(index >= 0 && index < Q) ArrBin[index] += 1;
                    }
                }
            }
        }
        for(int k = 0; k < Q; k++) {
            fprintf(filehisto, "%s %i %i\n", NomiMetodi[metodo], k, ArrBin[k]);
        }
        free(ArrBin);
        free(Temp);
        free(Position);
    }
    fclose(filegraph);
    fclose(filefit);
    fclose(filehisto);
    return 0;
}

// Gestisce la chiamata al generatore specifico e il passo
void SM(struct PQ *S, int *seed, int method)
{
    double val_norm = 0.0;
    // Selezione del generatore e normalizzazione in [0,1]
    if(method == 0) {
        *seed = LEcuyer1(*seed);
        val_norm = (double)((unsigned int)*seed) / 4294967296.0;
    } else if(method == 1) {
        *seed = LEcuyer2(*seed);
        val_norm = (double)(*seed) / 2147483647.0;
    } else if(method == 2) {
        *seed = MINSTAND(*seed);
        val_norm = (double)(*seed) / 2147483647.0;
    } else if(method == 3) {
        *seed = INFRAND(*seed);
        val_norm = (double)(*seed) / 2147483648.0;
    }
    
    S->q = val_norm;
    // Decisione del passo
    if(S->q >= 0.5) S->p = 1;
    else S->p = -1; 
}

// L'Ecuyer1 
int LEcuyer1(int In) {
    unsigned long long a = 1181783497276652981ULL;
    return (int)(a*In);
}

// L'Ecuyer 2
int LEcuyer2(int In) {
    long long a = 1385320207;
    long long m = (1LL << 31) - 1;
    return (int)((a*In)%m);
}

// MinStand
int MINSTAND(int In) {
    int a = 16807;
    long long m = (1LL << 31) - 1;
    return (int)(((long long)a * In) % m);
}

// Infamous Randu
int INFRAND(int In) {
    int a = (1 << 16) + 3;
    long long m = 1LL << 31;
    return (int)(((long long)a * In) % m);
}
