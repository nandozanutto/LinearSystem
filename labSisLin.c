#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "SistemasLineares.h"



int main ()
{
    real_t *x;
    real_t *res;
    real_t norma;
    int i=0;
    int iteracoes=0;
    double *tTotal;
    long int pos;
    char c;
    

    while(1){
        SistLinear_t *SL = lerSistLinear();
        SistLinear_t *newSL = alocaSistLinear(SL->n);
        x = (real_t *)calloc(SL->n, sizeof(real_t));
        tTotal = (double *)malloc(sizeof(double));
        res = (real_t *)calloc(SL->n, sizeof(real_t));
        
        //***salvando sistema original
        for(int i=0; i<SL->n; i++)
            for(int j=0; j<SL->n; j++)
                newSL->A[i][j] = SL->A[i][j];
        for(int i=0; i<SL->n; i++)
            newSL->b[i] = SL->b[i];
        newSL->n = SL->n;
        newSL->erro = SL->erro;
        //*****

        i++;
        eliminacaoGauss(SL, x, tTotal); 
        printf("***** Sistema %d --> n = %d, erro: %f\n", i, SL->n, SL->erro);
        printf("===> Eliminação Gauss: %lf ms\n --> X: ", tTotal[0]);
        prnVetor(x, SL->n);
        norma = normaL2Residuo(newSL, x, res);
        printf(" --> Norma L2 do residuo: %.7g\n", norma);

        if(norma > 5.0){
            iteracoes = refinamento(newSL, x, tTotal);
            printf("\n===> Refinamento: %lf ms --> %d iteracoes\n --> X: ", tTotal[0], iteracoes);
            prnVetor(x, SL->n);
            norma = normaL2Residuo(newSL, x, res);
            printf(" --> Norma L2 do residuo: %.7g\n", norma);
        }

        iteracoes = gaussJacobi(newSL, x, tTotal);
        printf("\n===> Jacobi: %lf ms --> %d iteracoes\n --> X: ", tTotal[0], iteracoes);
        prnVetor(x, SL->n);
        norma = normaL2Residuo(newSL, x, res);
        printf(" --> Norma L2 do residuo: %.7g\n", norma);

        if(norma > 5.0){
            iteracoes = refinamento(newSL, x, tTotal);
            printf("\n===> Refinamento: %lf ms --> %d iteracoes\n --> X: ", tTotal[0], iteracoes);
            prnVetor(x, SL->n);
            norma = normaL2Residuo(newSL, x, res);
            printf(" --> Norma L2 do residuo: %.7g\n", norma);
        }
        
        iteracoes = gaussSeidel(newSL, x, tTotal);
        printf("\n===> Gauss-Seidel: %lf ms --> %d iteracoes\n --> X: ", tTotal[0], iteracoes);
        prnVetor(x, SL->n);
        norma = normaL2Residuo(newSL, x, res);
        printf(" --> Norma L2 do residuo: %.7g\n", norma);

        if(norma > 5.0){
            iteracoes = refinamento(newSL, x, tTotal);
            printf("\n===> Refinamento: %lf ms --> %d iteracoes\n --> X: ", tTotal[0], iteracoes);
            prnVetor(x, SL->n);
            norma = normaL2Residuo(newSL, x, res);
            printf(" --> Norma L2 do residuo: %.7g\n", norma);
        }
        printf("\n");

        pos = ftell(stdin);
        for(c=fgetc(stdin); c == '\n' || c == ' '; c=fgetc(stdin))//searching for new linearSystems
            pos = ftell(stdin);//saving position of pointer
        if(c == EOF) break;
        else fseek (stdin, pos, SEEK_SET);//Number found!!
    
        liberaSistLinear(SL);
        liberaSistLinear(newSL);
    }
    
    
}

