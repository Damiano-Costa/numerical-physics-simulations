#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <time.h>

/*
 * Il codice utilizza il metodo di BoxMuller per generare una gaussiana. 
 * Richiede in input un numero grande corrispondente al numero di coppie 
 * rho1 rho2 che andrà a generare. 
 * L'output viene scritto su un file txt chiamato boxmul.txt. 
 */

void generator(double *BM, int Ncouples);
void MaxMin(double *BM, int Ncouples, double arr[2]);

int main(int argc, char **argv)
{
    // Controllo argomenti
    if(argc != 2)
    {
        printf("errore negli argomenti\n");
        return 0;
    }
    
    // Esecuzione programma
    if(argc == 2)
    {
        srand(time((NULL)));
        
        FILE *file = fopen("boxmul.txt", "w");
        double MaMi[2] = {0};
        
        // Acquisizione input
        int Ncouples = atoi(argv[1]);
        int Nbins = Ncouples / 10;
        
        // Allocazione array
        double *BM = calloc(Ncouples, sizeof(double));
        
        // Generazione numeri casuali
        generator(BM, Ncouples);
        
        // Ricerca massimi e minimi
        MaxMin(BM, Ncouples, MaMi);
        
        double dif = MaMi[1] - MaMi[0];
        // MaMi[0] è il minimo; normalizzo rispetto al minimo
        for(int i = 0; i < Ncouples; i++) 
            BM[i] -= MaMi[0];
            
        // Calcolo istogramma
        int *freq = calloc(Nbins, sizeof(int));
        for(int i = 0; i < Ncouples; i++)
        {
            int bin = (int)((BM[i] / dif) * Nbins);
            if(bin == Nbins) bin = Nbins - 1;
            freq[bin] += 1;
        }
        
        // Scrittura su file
        for(int i = 0; i < Nbins; i++)
        {
            fprintf(file, "%i\n", freq[i]);
        }
        
        fclose(file);
    }
}

void generator(double *BM, int Ncouples)
{
    for(int k = 0; k < Ncouples; k++)
    {
        double q = (double)rand() / RAND_MAX;
        double p = (double)rand() / RAND_MAX;
        // Formula Box-Muller vista in classe
        double u = sqrt((1 - log(q))) * sin(2 * M_PI * p);
        BM[k] = u;
    }
}

void MaxMin(double *BM, int Ncouples, double arr[2])
{
    double max = BM[0];
    double min = BM[0];
    
    // Trova il massimo
    for(int i = 1; i < Ncouples; i++)
    {
        if(BM[i] > max) 
            max = BM[i];
    }
    
    // Trova il minimo
    for(int i = 1; i < Ncouples; i++)
    {
        if(BM[i] < min) 
            min = BM[i];
    }
    
    arr[0] = min;
    arr[1] = max;
}
