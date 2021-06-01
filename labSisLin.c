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
    double *tTotal;
    long int pos;
    char c;
    

    while(1){
        SistLinear_t *SL = lerSistLinear();
        SistLinear_t *newSL = alocaSistLinear(SL->n);
        x = (real_t *)calloc(SL->n, sizeof(real_t));
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
        printf("===> Eliminação Gauss: tempo\n --> X: ");
        prnVetor(x, SL->n);
        norma = normaL2Residuo(newSL, x, res);
        printf(" --> Norma L2 do residuo: %.7g\n", norma);

        if(norma > 5.0){
            printf("\n===> Refinamento: tempo --> %d iteracoes\n --> X: ", refinamento(newSL, x, tTotal));
            prnVetor(x, SL->n);
            norma = normaL2Residuo(newSL, x, res);
            printf(" --> Norma L2 do residuo: %.7g\n", norma);
        }

        printf("\n===> Jacobi: tempo --> %d iteracoes\n --> X: ", gaussJacobi(newSL, x, tTotal));
        prnVetor(x, SL->n);
        norma = normaL2Residuo(newSL, x, res);
        printf(" --> Norma L2 do residuo: %.7g\n", norma);

        if(norma > 5.0){
            printf("\n===> Refinamento: tempo --> %d iteracoes\n --> X: ", refinamento(newSL, x, tTotal));
            prnVetor(x, SL->n);
            norma = normaL2Residuo(newSL, x, res);
            printf(" --> Norma L2 do residuo: %.7g\n", norma);
        }
        
        printf("\n===> Gauss-Seidel: tempo --> %d iteracoes\n --> X: ", gaussSeidel(newSL, x, tTotal));
        prnVetor(x, SL->n);
        norma = normaL2Residuo(newSL, x, res);
        printf(" --> Norma L2 do residuo: %.7g\n", norma);

        if(norma > 5.0){
            printf("\n===> Refinamento: tempo --> %d iteracoes\n --> X: ", refinamento(newSL, x, tTotal));
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
    
    }
    
    
}

